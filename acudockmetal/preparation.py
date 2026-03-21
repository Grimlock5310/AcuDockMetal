"""
Receptor and Ligand Preparation Pipeline.

Handles structure repair (PDBFixer), protonation, metal-site remediation,
box definition, ligand conformer generation, tautomer enumeration, and
format conversion for docking engines (PDBQT via Meeko).
"""

from __future__ import annotations

import logging
import os
import tempfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

import numpy as np

from acudockmetal.chelator_detect import ChelatingGroupDetector, DonorGroup
from acudockmetal.metal_site import MetalSite, MetalSiteDetector

log = logging.getLogger(__name__)

# Elements recognized as metal ions in PDB HETATM records
_METAL_ELEMENTS = {
    "ZN", "FE", "CU", "CO", "NI", "MN", "MG", "CA", "NA", "K",
    "LA", "CE", "PR", "ND", "SM", "EU", "GD", "TB",
    "DY", "HO", "ER", "TM", "YB", "LU",
}


@dataclass
class PreparedReceptor:
    """Container for a receptor prepared for docking."""
    pdb_path: str              # path to cleaned PDB
    pdbqt_path: Optional[str]  # path to PDBQT (for Vina)
    metal_sites: List[MetalSite]
    box_center: np.ndarray     # (3,)
    box_size: np.ndarray       # (3,)  half-widths in angstrom
    ph: float
    notes: List[str] = field(default_factory=list)


@dataclass
class PreparedLigand:
    """Container for a prepared ligand variant."""
    mol: object                # RDKit Mol
    smiles: str
    name: str
    pdbqt_string: Optional[str]
    pdb_path: Optional[str]
    donor_groups: List[DonorGroup]
    variant_label: str         # e.g. "tautomer_0", "protomer_1"

    @property
    def total_denticity(self) -> int:
        return sum(d.denticity for d in self.donor_groups)


class ReceptorPreparator:
    """
    Prepare a protein receptor for docking.

    Steps:
    1. Repair structure with PDBFixer (missing atoms, residues)
    2. Protonate at specified pH
    3. Detect and characterize metal sites
    4. Define docking box (centered on metal site or user-specified)
    5. Convert to PDBQT via Meeko
    """

    def __init__(
        self,
        ph: float = 7.4,
        box_padding: float = 10.0,
        keep_waters_near_metal: bool = True,
        water_distance_cutoff: float = 3.5,
        detect_metals: bool = True,
    ) -> None:
        self.ph = ph
        self.box_padding = box_padding
        self.keep_waters_near_metal = keep_waters_near_metal
        self.water_distance_cutoff = water_distance_cutoff
        self.detect_metals = detect_metals

    def prepare(
        self,
        pdb_path: str,
        box_center: Optional[np.ndarray] = None,
        box_size: Optional[np.ndarray] = None,
        output_dir: Optional[str] = None,
    ) -> PreparedReceptor:
        """
        Full receptor preparation pipeline.

        Args:
            pdb_path: input PDB file.
            box_center: optional explicit box center (auto-detected from
                        metal site if not given).
            box_size: optional explicit box half-widths.
            output_dir: directory for output files (temp dir if None).

        Returns:
            PreparedReceptor with paths, metal sites, and box info.
        """
        if output_dir is None:
            output_dir = tempfile.mkdtemp(prefix="acudock_")
        os.makedirs(output_dir, exist_ok=True)

        notes: List[str] = []

        # --- Step 1: Repair with PDBFixer ---
        fixed_path = os.path.join(output_dir, "receptor_fixed.pdb")
        try:
            self._fix_structure(pdb_path, fixed_path)
            notes.append("Structure repaired with PDBFixer.")
        except Exception as e:
            log.warning("PDBFixer failed (%s); using original structure.", e)
            fixed_path = pdb_path
            notes.append(f"PDBFixer skipped: {e}")

        # --- Step 2: Detect metal sites ---
        metal_sites: List[MetalSite] = []
        if self.detect_metals:
            detector = MetalSiteDetector()
            # Try detection on fixed structure first
            metal_sites = detector.detect_from_pdb(fixed_path)
            # Fallback: detect on ORIGINAL PDB (PDBFixer/OpenMM may strip
            # element columns or rename residues, breaking detection)
            if not metal_sites and fixed_path != pdb_path:
                metal_sites = detector.detect_from_pdb(pdb_path)
                if metal_sites:
                    notes.append(
                        "Metal sites detected from original PDB "
                        "(PDBFixer output lacked proper metal records).")
            notes.append(f"Detected {len(metal_sites)} metal site(s).")

        # --- Step 3: Auto-define box ---
        if box_center is None:
            if metal_sites:
                box_center = metal_sites[0].coord.copy()
                notes.append(
                    f"Box centered on metal site: {metal_sites[0].label}")
            else:
                # Fallback: center of mass of all atoms
                box_center = self._center_of_mass(fixed_path)
                notes.append("Box centered on protein center of mass.")

        if box_size is None:
            box_size = np.array([self.box_padding] * 3)

        # --- Step 4: Convert to PDBQT ---
        pdbqt_path = os.path.join(output_dir, "receptor.pdbqt")
        try:
            self._to_pdbqt(fixed_path, pdbqt_path)
            notes.append("PDBQT generated via Meeko/obabel.")
        except Exception as e:
            log.warning("PDBQT conversion failed: %s", e)
            pdbqt_path = None
            notes.append(f"PDBQT conversion skipped: {e}")

        return PreparedReceptor(
            pdb_path=fixed_path,
            pdbqt_path=pdbqt_path,
            metal_sites=metal_sites,
            box_center=box_center,
            box_size=box_size,
            ph=self.ph,
            notes=notes,
        )

    def _fix_structure(self, input_path: str, output_path: str) -> None:
        """Use PDBFixer to repair missing atoms/residues and protonate."""
        from pdbfixer import PDBFixer
        from openmm.app import PDBFile

        fixer = PDBFixer(filename=input_path)
        fixer.findMissingResidues()
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()
        fixer.addMissingHydrogens(self.ph)

        if not self.keep_waters_near_metal:
            fixer.removeHeterogens(keepWater=False)

        with open(output_path, "w") as f:
            PDBFile.writeFile(fixer.topology, fixer.positions, f)

        # PDBFixer may strip metal ions — re-inject them if missing
        self._restore_metals(input_path, output_path)

    def _restore_metals(self, original_pdb: str, fixed_pdb: str) -> None:
        """Re-inject metal HETATM lines if PDBFixer dropped them."""
        original_metals = []
        with open(original_pdb) as f:
            for line in f:
                if line.startswith("HETATM"):
                    elem = line[76:78].strip().upper() if len(line) >= 78 else ""
                    resname = line[17:20].strip().upper()
                    # Check both element column and residue name
                    if elem in _METAL_ELEMENTS or resname in _METAL_ELEMENTS:
                        original_metals.append(line)

        if not original_metals:
            log.info("No metal HETATM lines found in original PDB %s",
                     original_pdb)
            return

        log.info("Found %d metal HETATM line(s) in original PDB",
                 len(original_metals))

        with open(fixed_pdb) as f:
            fixed_lines = f.readlines()

        metals_present = any(
            (line[76:78].strip().upper() if len(line) >= 78 else "")
            in _METAL_ELEMENTS
            or line[17:20].strip().upper() in _METAL_ELEMENTS
            for line in fixed_lines
            if line.startswith("HETATM")
        )

        if metals_present:
            log.info("Metals already present in fixed PDB — no re-injection needed")
            return

        log.info("Re-injecting %d metal ion(s) stripped by PDBFixer",
                 len(original_metals))
        inserted = False
        with open(fixed_pdb, "w") as f:
            for line in fixed_lines:
                stripped = line.strip()
                if stripped in ("END", "ENDMDL"):
                    for m in original_metals:
                        f.write(m if m.endswith("\n") else m + "\n")
                    inserted = True
                f.write(line)
            if not inserted:
                # No END record found — append metals at end of file
                for m in original_metals:
                    f.write(m if m.endswith("\n") else m + "\n")

    def _to_pdbqt(self, pdb_path: str, pdbqt_path: str) -> None:
        """Convert PDB to PDBQT for AutoDock Vina."""
        import shutil
        import subprocess

        # Primary: mk_prepare_receptor.py (installed with meeko)
        mk_prep = shutil.which("mk_prepare_receptor.py")
        if mk_prep:
            out_dir = os.path.dirname(pdbqt_path) or "."
            basename = os.path.splitext(os.path.basename(pdbqt_path))[0]
            abs_pdb = os.path.abspath(pdb_path)
            result = subprocess.run(
                [mk_prep, "--read_pdb", abs_pdb,
                 "-o", basename, "-p", "--allow_bad_res"],
                capture_output=True, text=True, timeout=300,
                cwd=out_dir,
            )
            if result.returncode == 0:
                expected = os.path.join(out_dir, f"{basename}.pdbqt")
                if os.path.exists(expected) and expected != pdbqt_path:
                    shutil.move(expected, pdbqt_path)
                return
            log.warning("mk_prepare_receptor failed: %s",
                        result.stderr.strip()[:500])

        # Fallback: obabel
        result = subprocess.run(
            ["obabel", pdb_path, "-O", pdbqt_path,
             "-xr", "--partialcharge", "gasteiger"],
            capture_output=True, text=True, timeout=120,
        )
        if result.returncode != 0:
            raise RuntimeError(
                "PDBQT conversion failed — neither "
                "mk_prepare_receptor.py nor obabel available")

    def _center_of_mass(self, pdb_path: str) -> np.ndarray:
        """Compute center of mass from a PDB file."""
        coords = []
        with open(pdb_path) as f:
            for line in f:
                if line.startswith(("ATOM", "HETATM")):
                    try:
                        x = float(line[30:38])
                        y = float(line[38:46])
                        z = float(line[46:54])
                        coords.append([x, y, z])
                    except ValueError:
                        continue
        if not coords:
            return np.zeros(3)
        return np.mean(coords, axis=0)


class LigandPreparator:
    """
    Prepare ligands for metal-aware docking.

    Steps:
    1. Parse SMILES or SDF; sanitize with RDKit
    2. Add hydrogens and generate 3D conformers
    3. Detect chelating groups
    4. Enumerate tautomers/protomers relevant to metal binding
    5. Convert to PDBQT via Meeko
    """

    def __init__(
        self,
        num_conformers: int = 10,
        enumerate_tautomers: bool = True,
        max_tautomers: int = 5,
        chelator_min_confidence: float = 0.3,
    ) -> None:
        self.num_conformers = num_conformers
        self.enumerate_tautomers = enumerate_tautomers
        self.max_tautomers = max_tautomers
        self.chelator_detector = ChelatingGroupDetector(
            min_confidence=chelator_min_confidence)

    def prepare_from_smiles(
        self,
        smiles: str,
        name: str = "ligand",
        output_dir: Optional[str] = None,
    ) -> List[PreparedLigand]:
        """Prepare a ligand from SMILES string."""
        from rdkit import Chem
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError(f"Invalid SMILES: {smiles}")
        mol.SetProp("_Name", name)
        return self._prepare_mol(mol, name, output_dir, original_smiles=smiles)

    def prepare_from_sdf(
        self,
        sdf_path: str,
        output_dir: Optional[str] = None,
    ) -> List[PreparedLigand]:
        """Prepare ligands from an SDF file (first molecule)."""
        from rdkit import Chem
        suppl = Chem.SDMolSupplier(sdf_path, removeHs=False)
        results = []
        for i, mol in enumerate(suppl):
            if mol is None:
                continue
            name = mol.GetProp("_Name") if mol.HasProp("_Name") else f"lig_{i}"
            results.extend(self._prepare_mol(mol, name, output_dir))
        return results

    def _prepare_mol(
        self,
        mol,
        name: str,
        output_dir: Optional[str],
        original_smiles: Optional[str] = None,
    ) -> List[PreparedLigand]:
        """Core preparation for a single molecule."""
        from rdkit import Chem
        from rdkit.Chem import AllChem

        if output_dir is None:
            output_dir = tempfile.mkdtemp(prefix="acudock_lig_")
        os.makedirs(output_dir, exist_ok=True)

        mol = Chem.AddHs(mol)

        # Check for metal atoms (affects conformer gen + PDBQT path)
        has_metal = any(
            a.GetAtomicNum() > 20
            or a.GetAtomicNum() in (3, 4, 11, 12, 13, 19, 20)
            for a in mol.GetAtoms()
        )

        # Generate 3D conformer
        params = AllChem.ETKDGv3()
        params.numThreads = 0
        params.randomSeed = 42
        conf_ok = AllChem.EmbedMolecule(mol, params)
        if conf_ok == -1:
            conf_ok = AllChem.EmbedMolecule(mol, AllChem.ETKDG())
        if conf_ok == -1 and has_metal:
            # Metals: distance geometry with random coords (no FF needed)
            params2 = AllChem.ETKDGv3()
            params2.useRandomCoords = True
            params2.randomSeed = 42
            conf_ok = AllChem.EmbedMolecule(mol, params2)
        if conf_ok != -1:
            try:
                AllChem.MMFFOptimizeMolecule(mol, maxIters=500)
            except Exception:
                pass  # MMFF may not support metals
        else:
            log.warning("Could not generate 3D conformer; "
                        "PDBQT conversion may fail.")

        variants = []

        # --- Tautomer enumeration ---
        mols_to_process = [("base", mol)]
        if self.enumerate_tautomers:
            try:
                from rdkit.Chem.MolStandardize import rdMolStandardize
                enumerator = rdMolStandardize.TautomerEnumerator()
                tauts = enumerator.Enumerate(mol)
                for ti, t in enumerate(tauts):
                    if ti >= self.max_tautomers:
                        break
                    t = Chem.AddHs(t)
                    AllChem.EmbedMolecule(t, AllChem.ETKDGv3())
                    AllChem.MMFFOptimizeMolecule(t, maxIters=200)
                    mols_to_process.append((f"tautomer_{ti}", t))
            except Exception as e:
                log.warning("Tautomer enumeration failed: %s", e)

        # --- Process each variant ---
        for label, m in mols_to_process:
            donors = self.chelator_detector.detect(m)
            smi = Chem.MolToSmiles(Chem.RemoveHs(m))

            # Write PDB
            pdb_path = os.path.join(output_dir, f"{name}_{label}.pdb")
            Chem.MolToPDBFile(m, pdb_path)

            # Convert to PDBQT
            pdbqt_str = None
            if has_metal:
                # Meeko can't handle metals — use obabel directly
                # Skip Gasteiger charges: they silently drop metal complexes
                try:
                    pdbqt_str = self._to_pdbqt_obabel(pdb_path,
                                                       use_gasteiger=False)
                except Exception as e:
                    log.warning("obabel PDBQT failed for metal ligand %s: %s",
                                label, e)
                # If RDKit conformer failed (degenerate PDB), generate 3D
                # coords and PDBQT directly from SMILES via obabel --gen3d
                if not pdbqt_str or not pdbqt_str.strip():
                    # Use original SMILES (not RDKit canonical) — RDKit may
                    # mangle metal coordination in canonical form
                    gen3d_smi = original_smiles or smi
                    try:
                        pdbqt_str = self._smiles_to_pdbqt_obabel(gen3d_smi)
                    except Exception as e2:
                        log.warning(
                            "obabel --gen3d from SMILES also failed for %s: %s",
                            label, e2)
            else:
                try:
                    pdbqt_str = self._to_pdbqt_meeko(m)
                except Exception as e:
                    log.warning("Meeko PDBQT conversion failed for %s: %s",
                                label, e)
                    # Fallback: try obabel
                    try:
                        pdbqt_str = self._to_pdbqt_obabel(pdb_path)
                    except Exception as e2:
                        log.warning("obabel PDBQT fallback failed: %s", e2)

            variants.append(PreparedLigand(
                mol=m,
                smiles=smi,
                name=name,
                pdbqt_string=pdbqt_str,
                pdb_path=pdb_path,
                donor_groups=donors,
                variant_label=label,
            ))

        return variants

    def _to_pdbqt_meeko(self, mol) -> str:
        """Convert RDKit Mol to PDBQT string using Meeko."""
        from meeko import MoleculePreparation, PDBQTWriterLegacy

        preparator = MoleculePreparation()
        mol_setup_list = preparator.prepare(mol)
        if not mol_setup_list:
            raise RuntimeError("Meeko preparation returned empty result")
        pdbqt_str, is_ok, err = PDBQTWriterLegacy.write_string(
            mol_setup_list[0])
        if not is_ok:
            raise RuntimeError(f"Meeko PDBQT write error: {err}")
        return pdbqt_str

    def _to_pdbqt_obabel(self, pdb_path: str,
                          use_gasteiger: bool = True) -> str:
        """Fallback PDBQT conversion using Open Babel."""
        import subprocess
        cmd = ["obabel", pdb_path, "-O", "-", "-opdbqt"]
        # Gasteiger charges silently drop metal complexes in obabel 3.x
        if use_gasteiger:
            cmd += ["--partialcharge", "gasteiger"]
        result = subprocess.run(
            cmd, capture_output=True, text=True, timeout=60,
        )
        if result.returncode != 0:
            raise RuntimeError(f"obabel failed: {result.stderr.strip()}")
        if not result.stdout.strip():
            raise RuntimeError("obabel returned empty PDBQT")
        return result.stdout

    def _smiles_to_pdbqt_obabel(self, smiles: str) -> str:
        """Generate 3D coords and PDBQT directly from SMILES via obabel.

        This bypasses RDKit conformer generation entirely, which is
        necessary for metal-containing ligands where RDKit's distance
        geometry fails. Does NOT use Gasteiger charges — they silently
        drop metal complexes in Open Babel 3.x.
        """
        import subprocess
        result = subprocess.run(
            ["obabel", "-ismi", "-opdbqt", "--gen3d"],
            input=smiles,
            capture_output=True, text=True, timeout=120,
        )
        if result.returncode != 0:
            raise RuntimeError(
                f"obabel --gen3d failed: {result.stderr.strip()}")
        if not result.stdout.strip():
            raise RuntimeError("obabel --gen3d returned empty PDBQT")
        return result.stdout

    def summary(self, variants: List[PreparedLigand]) -> str:
        """Human-readable summary of prepared ligand variants."""
        lines = [f"Prepared {len(variants)} ligand variant(s):"]
        for v in variants:
            n_donors = len(v.donor_groups)
            dent = v.total_denticity
            has_pdbqt = "yes" if v.pdbqt_string else "no"
            lines.append(
                f"  {v.name}/{v.variant_label}: {v.smiles[:50]} "
                f"donors={n_donors} denticity={dent} pdbqt={has_pdbqt}"
            )
        return "\n".join(lines)
