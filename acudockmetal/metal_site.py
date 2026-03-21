"""
Metal Site Detection and Coordination Hypothesis Generation.

Detects metal ions in protein structures (PDB/mmCIF), identifies
coordinating residues, classifies coordination geometry, and enumerates
plausible coordination hypotheses for docking.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import numpy as np

from acudockmetal.metal_params import (
    GEOMETRY_TEMPLATES,
    GeometryTemplate,
    MetalParameterLibrary,
)

# Known metal ion residue names in PDB files
METAL_RESIDUES = {
    "ZN", "FE", "FE2", "CU", "CU1", "CO", "NI", "MN",
    "MG", "CA", "NA", "K",
    "LA", "CE", "PR", "ND", "SM", "EU", "GD", "TB",
    "DY", "HO", "ER", "TM", "YB", "LU",
}

# Element symbols for fallback detection (handles PDBFixer renaming)
_METAL_SYMBOLS = {
    "ZN", "FE", "CU", "CO", "NI", "MN", "MG", "CA", "NA", "K",
    "LA", "CE", "PR", "ND", "SM", "EU", "GD", "TB",
    "DY", "HO", "ER", "TM", "YB", "LU",
}

# Map PDB residue names to our library keys
_RESNAME_TO_SYMBOL = {
    "ZN": "ZN", "FE": "FE", "FE2": "FE", "CU": "CU", "CU1": "CU",
    "CO": "CO", "NI": "NI", "MN": "MN", "MG": "MG", "CA": "CA",
    # Lanthanides all map to generic "LN"
    "LA": "LN", "CE": "LN", "PR": "LN", "ND": "LN", "SM": "LN",
    "EU": "LN", "GD": "LN", "TB": "LN", "DY": "LN", "HO": "LN",
    "ER": "LN", "TM": "LN", "YB": "LN", "LU": "LN",
}


@dataclass
class AtomInfo:
    """Lightweight descriptor for a coordinating atom."""
    chain_id: str
    residue_name: str
    residue_number: int
    atom_name: str
    element: str
    coord: np.ndarray  # (3,)
    is_water: bool = False

    def label(self) -> str:
        chain = self.chain_id.strip() or "?"
        return f"{chain}/{self.residue_name}{self.residue_number}/{self.atom_name}"


@dataclass
class MetalSite:
    """A detected metal-ion site in a protein structure."""
    metal_symbol: str          # library key (ZN, FE, …)
    pdb_resname: str           # original PDB residue name
    chain_id: str
    residue_number: int
    coord: np.ndarray          # (3,) metal position
    coordinating_atoms: List[AtomInfo] = field(default_factory=list)
    observed_cn: int = 0
    classified_geometry: Optional[str] = None
    geometry_rmsd: float = float("inf")
    has_incomplete_coordination: bool = False

    @property
    def label(self) -> str:
        chain = self.chain_id.strip() or "?"
        return f"{self.metal_symbol}({chain}/{self.pdb_resname}{self.residue_number})"


@dataclass
class CoordinationHypothesis:
    """
    A discrete hypothesis about the coordination state of a metal site.

    The docking system docks the ligand independently under each hypothesis
    and compares results across hypotheses.
    """
    hypothesis_id: int
    metal_site: MetalSite
    metal_symbol: str
    oxidation_state: int
    coordination_number: int
    geometry: str              # geometry template name
    protein_donors: List[AtomInfo]
    water_positions: List[np.ndarray]
    open_slots: int            # slots available for ligand binding
    weight: float              # prior probability / ranking weight

    @property
    def label(self) -> str:
        return (f"H{self.hypothesis_id}: "
                f"{self.metal_symbol}(+{self.oxidation_state}) "
                f"CN={self.coordination_number} "
                f"{self.geometry} "
                f"open={self.open_slots}")


class MetalSiteDetector:
    """
    Detect and characterize metal-ion sites in a protein structure.

    Reads a PDB file using BioPython, locates metal ions, identifies
    coordinating atoms within distance cutoffs, and classifies the
    observed coordination geometry.
    """

    def __init__(
        self,
        distance_cutoff: float = 2.8,
        metal_library: Optional[MetalParameterLibrary] = None,
    ) -> None:
        self.distance_cutoff = distance_cutoff
        self.metal_lib = metal_library or MetalParameterLibrary()

    def detect_from_pdb(self, pdb_path: str) -> List[MetalSite]:
        """
        Detect metal sites in a PDB file.

        Args:
            pdb_path: path to a PDB-format file.

        Returns:
            List of MetalSite objects.
        """
        from Bio.PDB import PDBParser

        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("protein", pdb_path)

        # Collect all atoms for neighbor search
        all_atoms = []
        all_coords = []
        atom_info_map = {}
        metal_atoms = []

        for model in structure:
            for chain in model:
                for residue in chain:
                    resname = residue.get_resname().strip()
                    resnum = residue.get_id()[1]
                    chain_id = chain.get_id()
                    is_water = resname in ("HOH", "WAT", "H2O", "DOD")

                    for atom in residue:
                        coord = atom.get_vector().get_array()
                        elem = atom.element.strip().upper() if atom.element else ""
                        # Fallback: infer element from atom name when
                        # element column is missing (common with OpenMM output)
                        if not elem:
                            aname = atom.get_name().strip().upper()
                            if aname in _METAL_SYMBOLS:
                                elem = aname
                        ainfo = AtomInfo(
                            chain_id=chain_id,
                            residue_name=resname,
                            residue_number=resnum,
                            atom_name=atom.get_name(),
                            element=elem,
                            coord=np.array(coord, dtype=np.float64),
                            is_water=is_water,
                        )
                        idx = len(all_atoms)
                        all_atoms.append(ainfo)
                        all_coords.append(coord)

                        if resname in METAL_RESIDUES or elem in _METAL_SYMBOLS:
                            metal_atoms.append(idx)

        if not all_coords:
            return []

        all_coords_arr = np.array(all_coords)
        sites: List[MetalSite] = []

        for mi in metal_atoms:
            metal_info = all_atoms[mi]
            metal_coord = all_coords_arr[mi]
            metal_sym = _RESNAME_TO_SYMBOL.get(
                metal_info.residue_name, metal_info.element)

            # Find coordinating atoms within cutoff
            dists = np.linalg.norm(all_coords_arr - metal_coord, axis=1)
            nearby_mask = (dists < self.distance_cutoff) & (dists > 0.5)
            nearby_indices = np.where(nearby_mask)[0]

            coord_atoms = []
            for ni in nearby_indices:
                ai = all_atoms[ni]
                # Only consider N, O, S as potential donors
                if ai.element in ("N", "O", "S"):
                    coord_atoms.append(ai)

            # Classify geometry
            geom_name, geom_rmsd = self._classify_geometry(
                metal_coord, coord_atoms)

            incomplete = False
            if metal_sym in self.metal_lib:
                mp = self.metal_lib[metal_sym]
                cn_lo, cn_hi = mp.typical_cn_range
                if len(coord_atoms) < cn_lo:
                    incomplete = True

            site = MetalSite(
                metal_symbol=metal_sym,
                pdb_resname=metal_info.residue_name,
                chain_id=metal_info.chain_id,
                residue_number=metal_info.residue_number,
                coord=metal_coord,
                coordinating_atoms=coord_atoms,
                observed_cn=len(coord_atoms),
                classified_geometry=geom_name,
                geometry_rmsd=geom_rmsd,
                has_incomplete_coordination=incomplete,
            )
            sites.append(site)

        return sites

    def _classify_geometry(
        self,
        metal_coord: np.ndarray,
        coord_atoms: List[AtomInfo],
    ) -> Tuple[Optional[str], float]:
        """
        Classify the observed coordination geometry by comparing
        observed angles to ideal geometry templates (FindGeo-style).
        """
        n = len(coord_atoms)
        if n < 2:
            return None, float("inf")

        # Compute all pairwise donor-metal-donor angles
        observed_angles = []
        for i in range(n):
            for j in range(i + 1, n):
                v1 = coord_atoms[i].coord - metal_coord
                v2 = coord_atoms[j].coord - metal_coord
                cos_a = np.dot(v1, v2) / (
                    np.linalg.norm(v1) * np.linalg.norm(v2) + 1e-12)
                cos_a = np.clip(cos_a, -1.0, 1.0)
                angle = np.degrees(np.arccos(cos_a))
                observed_angles.append(angle)

        # Compare against templates with matching CN
        best_name = None
        best_rmsd = float("inf")

        for name, tmpl in GEOMETRY_TEMPLATES.items():
            if tmpl.coordination_number != n:
                continue
            rmsd = tmpl.angle_deviation(observed_angles)
            if rmsd < best_rmsd:
                best_rmsd = rmsd
                best_name = name

        # Also try templates with CN +/- 1 for incomplete sites
        for name, tmpl in GEOMETRY_TEMPLATES.items():
            if abs(tmpl.coordination_number - n) == 1:
                rmsd = tmpl.angle_deviation(observed_angles)
                if rmsd < best_rmsd * 0.7:  # significantly better
                    best_rmsd = rmsd
                    best_name = name

        return best_name, best_rmsd

    def summary(self, sites: List[MetalSite]) -> str:
        """Human-readable summary of detected metal sites."""
        if not sites:
            return "No metal sites detected."
        lines = [f"Detected {len(sites)} metal site(s):"]
        for s in sites:
            donors = ", ".join(a.label() for a in s.coordinating_atoms
                               if not a.is_water)
            waters = sum(1 for a in s.coordinating_atoms if a.is_water)
            incomp = " [INCOMPLETE]" if s.has_incomplete_coordination else ""
            lines.append(
                f"  {s.label}: CN={s.observed_cn}, "
                f"geometry={s.classified_geometry or '?'} "
                f"(RMSD={s.geometry_rmsd:.1f} deg), "
                f"waters={waters}, donors=[{donors}]{incomp}"
            )
        return "\n".join(lines)


class CoordinationHypothesisGenerator:
    """
    Generate discrete coordination hypotheses for a metal site.

    Each hypothesis represents a plausible coordination state (CN, geometry,
    water occupancy, open slots for ligand) and is weighted by structural
    statistics.
    """

    def __init__(
        self,
        metal_library: Optional[MetalParameterLibrary] = None,
        max_hypotheses_per_site: int = 12,
    ) -> None:
        self.metal_lib = metal_library or MetalParameterLibrary()
        self.max_hypotheses = max_hypotheses_per_site

    def generate(self, site: MetalSite) -> List[CoordinationHypothesis]:
        """
        Generate coordination hypotheses for a metal site.

        Enumerates CN variants, geometry variants, and water occupancy
        states; weights by structural plausibility.
        """
        metal_sym = site.metal_symbol
        if metal_sym not in self.metal_lib:
            warnings.warn(f"Metal {metal_sym} not in library; "
                          f"generating minimal hypotheses.")
            return self._minimal_hypotheses(site)

        mp = self.metal_lib[metal_sym]
        protein_donors = [a for a in site.coordinating_atoms
                          if not a.is_water]
        n_protein = len(protein_donors)
        water_atoms = [a for a in site.coordinating_atoms if a.is_water]

        hypotheses: List[CoordinationHypothesis] = []
        hyp_id = 0

        for ox_state in mp.common_oxidation_states:
            for cn, geom_names in mp.preferred_geometries.items():
                for geom_name in geom_names:
                    # How many slots does the ligand get?
                    # open_slots = CN - protein_donors - waters_retained
                    max_open = cn - n_protein
                    if max_open < 1:
                        continue  # no room for ligand

                    # Enumerate water retention: 0 to min(max_open-1, n_waters)
                    max_waters = min(max_open - 1, len(water_atoms),
                                     cn - n_protein - 1)
                    max_waters = max(max_waters, 0)

                    for n_waters in range(max_waters + 1):
                        open_slots = cn - n_protein - n_waters
                        if open_slots < 1:
                            continue

                        # Compute weight (prior probability)
                        weight = self._compute_weight(
                            site, mp, cn, geom_name, n_waters, open_slots)

                        water_pos = [w.coord for w in water_atoms[:n_waters]]

                        hypotheses.append(CoordinationHypothesis(
                            hypothesis_id=hyp_id,
                            metal_site=site,
                            metal_symbol=metal_sym,
                            oxidation_state=ox_state,
                            coordination_number=cn,
                            geometry=geom_name,
                            protein_donors=protein_donors,
                            water_positions=water_pos,
                            open_slots=open_slots,
                            weight=weight,
                        ))
                        hyp_id += 1

        # Sort by weight and truncate
        hypotheses.sort(key=lambda h: h.weight, reverse=True)
        hypotheses = hypotheses[: self.max_hypotheses]

        # Renumber
        for i, h in enumerate(hypotheses):
            h.hypothesis_id = i

        return hypotheses

    def _compute_weight(
        self,
        site: MetalSite,
        mp,
        cn: int,
        geometry: str,
        n_waters: int,
        open_slots: int,
    ) -> float:
        """
        Compute a prior weight for a hypothesis based on:
        - match to observed CN
        - match to classified geometry
        - water propensity of the metal
        - structural plausibility
        """
        w = 1.0

        # Reward matching observed CN
        if cn == site.observed_cn:
            w *= 2.0
        elif abs(cn - site.observed_cn) == 1:
            w *= 1.0
        else:
            w *= 0.3

        # Reward matching classified geometry
        if geometry == site.classified_geometry:
            w *= 2.0

        # Water propensity
        if n_waters > 0:
            w *= mp.water_propensity ** n_waters
        else:
            w *= (1.0 - mp.water_propensity * 0.5)

        # Prefer fewer open slots (more constrained = more likely correct)
        w *= 1.0 / (1.0 + 0.2 * open_slots)

        # Default CN bonus
        if cn == mp.default_cn:
            w *= 1.5

        return w

    def _minimal_hypotheses(self, site: MetalSite) -> List[CoordinationHypothesis]:
        """Fallback hypotheses when metal is not in library."""
        protein_donors = [a for a in site.coordinating_atoms
                          if not a.is_water]
        cn = max(site.observed_cn, len(protein_donors) + 1)

        # Guess geometry from CN
        geom_map = {2: "linear", 3: "trigonal_planar", 4: "tetrahedral",
                    5: "trigonal_bipyramidal", 6: "octahedral"}
        geom = geom_map.get(cn, "octahedral")

        return [CoordinationHypothesis(
            hypothesis_id=0,
            metal_site=site,
            metal_symbol=site.metal_symbol,
            oxidation_state=2,
            coordination_number=cn,
            geometry=geom,
            protein_donors=protein_donors,
            water_positions=[],
            open_slots=max(cn - len(protein_donors), 1),
            weight=1.0,
        )]

    def summary(self, hypotheses: List[CoordinationHypothesis]) -> str:
        """Human-readable summary of hypotheses."""
        if not hypotheses:
            return "No hypotheses generated."
        lines = [f"Generated {len(hypotheses)} hypothesis/hypotheses:"]
        for h in hypotheses:
            lines.append(
                f"  {h.label} "
                f"(waters={len(h.water_positions)}, "
                f"weight={h.weight:.3f})"
            )
        return "\n".join(lines)
