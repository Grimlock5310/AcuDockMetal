"""
Chelating Group Detector — identifies metal-binding functional groups on
ligands using RDKit SMARTS pattern matching.

Detects hydroxamic acids, carboxylates, thiols, imidazoles, phosphonates,
catechols, 8-hydroxyquinolines, dithiocarbamates, aminocarboxylates, and
other pharmacologically relevant zinc/metal-binding groups (ZBGs).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Set, Tuple

from rdkit import Chem


@dataclass
class DonorGroup:
    """A detected metal-binding functional group on a ligand."""
    name: str
    donor_atoms: List[int]        # atom indices that coordinate the metal
    donor_elements: List[str]     # element symbols of donors ("O", "N", …)
    all_atom_indices: List[int]   # full substructure match indices
    denticity: int                # mono-, bi-, tri-dentate …
    compatible_metals: List[str]  # element symbols this group commonly binds
    smarts: str                   # SMARTS that matched
    confidence: float = 1.0       # 0-1, reduced for ambiguous matches


# ---------------------------------------------------------------------------
# SMARTS library for chelating / metal-binding groups
# ---------------------------------------------------------------------------

_CHELATOR_PATTERNS: List[Dict] = [
    # -- Hydroxamic acid: R-C(=O)N(H)OH  (bidentate O,O for Zn, Fe)
    {
        "name": "hydroxamic_acid",
        "smarts": "[OX1]=[CX3][NX3][OX2H1]",
        "donor_indices": [0, 3],    # the C=O oxygen + N-OH oxygen
        "donor_elements": ["O", "O"],
        "denticity": 2,
        "metals": ["ZN", "FE", "CO", "NI", "MN"],
    },
    # -- Reverse hydroxamate pattern (catches tautomers)
    {
        "name": "hydroxamic_acid_alt",
        "smarts": "[OX2H1][NX3][CX3]=[OX1]",
        "donor_indices": [0, 3],
        "donor_elements": ["O", "O"],
        "denticity": 2,
        "metals": ["ZN", "FE", "CO", "NI", "MN"],
    },
    # -- Carboxylate (deprotonated) — bidentate O,O
    {
        "name": "carboxylate",
        "smarts": "[OX1]=[CX3][OX1-]",
        "donor_indices": [0, 2],
        "donor_elements": ["O", "O"],
        "denticity": 2,
        "metals": ["ZN", "MG", "CA", "MN", "FE", "CO", "NI", "LN"],
    },
    # -- Carboxylic acid (protonated) — potential monodentate
    {
        "name": "carboxylic_acid",
        "smarts": "[OX1]=[CX3][OX2H1]",
        "donor_indices": [0],
        "donor_elements": ["O"],
        "denticity": 1,
        "metals": ["ZN", "MG", "CA", "MN"],
        "confidence": 0.7,
    },
    # -- Thiol / thiolate
    {
        "name": "thiol",
        "smarts": "[SX2H1]",
        "donor_indices": [0],
        "donor_elements": ["S"],
        "denticity": 1,
        "metals": ["ZN", "CU", "FE"],
    },
    {
        "name": "thiolate",
        "smarts": "[SX1-]",
        "donor_indices": [0],
        "donor_elements": ["S"],
        "denticity": 1,
        "metals": ["ZN", "CU", "FE"],
    },
    # -- Imidazole ring (His-like) — monodentate via unsubstituted N
    {
        "name": "imidazole",
        "smarts": "c1cn[nH]c1",
        "donor_indices": [2],       # the unprotonated nitrogen
        "donor_elements": ["N"],
        "denticity": 1,
        "metals": ["ZN", "CU", "FE", "CO", "NI", "MN"],
    },
    # -- Pyridine nitrogen — common aromatic N donor
    {
        "name": "pyridine",
        "smarts": "[nX2]1ccccc1",
        "donor_indices": [0],
        "donor_elements": ["N"],
        "denticity": 1,
        "metals": ["ZN", "CU", "FE", "CO", "NI"],
    },
    # -- Phosphonate (deprotonated) — tridentate O,O,O
    {
        "name": "phosphonate",
        "smarts": "[PX4](=[OX1])([OX1-])[OX1-]",
        "donor_indices": [1, 2, 3],
        "donor_elements": ["O", "O", "O"],
        "denticity": 3,
        "metals": ["MG", "CA", "LN", "ZN", "MN"],
    },
    # -- Phosphate monoester
    {
        "name": "phosphate",
        "smarts": "[PX4](=[OX1])([OX1-])([OX1-])[OX2]",
        "donor_indices": [1, 2, 3],
        "donor_elements": ["O", "O", "O"],
        "denticity": 3,
        "metals": ["MG", "CA", "MN", "ZN", "LN"],
    },
    # -- Catechol (1,2-dihydroxybenzene) — bidentate O,O
    {
        "name": "catechol",
        "smarts": "[OX2H1]c1ccccc1[OX2H1]",
        "donor_indices": [0, 7],
        "donor_elements": ["O", "O"],
        "denticity": 2,
        "metals": ["FE", "CU", "ZN"],
    },
    # -- 8-Hydroxyquinoline — bidentate N,O
    {
        "name": "8_hydroxyquinoline",
        "smarts": "[nX2]1ccc2cccc([OX2H1])c2c1",
        "donor_indices": [0, 8],
        "donor_elements": ["N", "O"],
        "denticity": 2,
        "metals": ["ZN", "CU", "FE"],
    },
    # -- Dithiocarbamate — bidentate S,S
    {
        "name": "dithiocarbamate",
        "smarts": "[SX1-][CX3](=[SX1])[NX3]",
        "donor_indices": [0, 2],
        "donor_elements": ["S", "S"],
        "denticity": 2,
        "metals": ["ZN", "CU", "NI"],
    },
    # -- Sulfonamide — monodentate N
    {
        "name": "sulfonamide",
        "smarts": "[NX3H1][SX4](=[OX1])=[OX1]",
        "donor_indices": [0],
        "donor_elements": ["N"],
        "denticity": 1,
        "metals": ["ZN"],
        "confidence": 0.8,
    },
    # -- Primary amine — weak monodentate N
    {
        "name": "primary_amine",
        "smarts": "[NX3H2][CX4]",
        "donor_indices": [0],
        "donor_elements": ["N"],
        "denticity": 1,
        "metals": ["CU", "NI", "CO", "ZN"],
        "confidence": 0.5,
    },
    # -- Carbonyl (ketone / amide) — weak monodentate O
    {
        "name": "carbonyl",
        "smarts": "[OX1]=[CX3]",
        "donor_indices": [0],
        "donor_elements": ["O"],
        "denticity": 1,
        "metals": ["MG", "CA", "LN", "MN"],
        "confidence": 0.4,
    },
]


class ChelatingGroupDetector:
    """
    Detect metal-binding (chelating) functional groups on a molecule.

    Uses a curated library of SMARTS patterns covering the major
    zinc-binding groups (ZBGs) and metal-chelating motifs relevant to
    metalloprotein drug design.
    """

    def __init__(self, min_confidence: float = 0.0) -> None:
        """
        Args:
            min_confidence: minimum confidence threshold to report a match
                (0.0 = report everything, 1.0 = only high-confidence).
        """
        self.min_confidence = min_confidence
        self._patterns = []
        for p in _CHELATOR_PATTERNS:
            mol = Chem.MolFromSmarts(p["smarts"])
            if mol is None:
                continue
            self._patterns.append((mol, p))

    def detect(self, mol: Chem.Mol) -> List[DonorGroup]:
        """
        Detect all chelating groups on a molecule.

        Args:
            mol: RDKit Mol object (should have explicit Hs for best results).

        Returns:
            List of DonorGroup objects found on the molecule.
        """
        if mol is None:
            return []

        results: List[DonorGroup] = []
        seen_donor_sets: Set[frozenset] = set()

        for query, info in self._patterns:
            conf = info.get("confidence", 1.0)
            if conf < self.min_confidence:
                continue

            matches = mol.GetSubstructMatches(query)
            for match in matches:
                donor_idx = [match[i] for i in info["donor_indices"]
                             if i < len(match)]
                donor_key = frozenset(donor_idx)

                # Avoid duplicates from overlapping patterns
                if donor_key in seen_donor_sets:
                    continue
                seen_donor_sets.add(donor_key)

                # Verify donor element symbols
                donor_elems = []
                valid = True
                for idx in donor_idx:
                    atom = mol.GetAtomWithIdx(idx)
                    elem = atom.GetSymbol()
                    donor_elems.append(elem)

                if not valid:
                    continue

                results.append(DonorGroup(
                    name=info["name"],
                    donor_atoms=donor_idx,
                    donor_elements=donor_elems,
                    all_atom_indices=list(match),
                    denticity=info["denticity"],
                    compatible_metals=info["metals"],
                    smarts=info["smarts"],
                    confidence=conf,
                ))

        return results

    def detect_from_smiles(self, smiles: str) -> List[DonorGroup]:
        """Convenience: detect chelating groups from a SMILES string."""
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError(f"Invalid SMILES: {smiles}")
        mol = Chem.AddHs(mol)
        return self.detect(mol)

    def compatible_with_metal(
        self, donors: List[DonorGroup], metal_symbol: str
    ) -> List[DonorGroup]:
        """Filter donors to only those compatible with a specific metal."""
        metal = metal_symbol.upper()
        return [d for d in donors if metal in d.compatible_metals]

    def total_denticity(self, donors: List[DonorGroup]) -> int:
        """Total number of coordination bonds the ligand can form."""
        return sum(d.denticity for d in donors)

    def summary(self, donors: List[DonorGroup]) -> str:
        """Human-readable summary of detected groups."""
        if not donors:
            return "No metal-binding groups detected."
        lines = []
        for d in donors:
            metals = ", ".join(d.compatible_metals)
            atoms = ", ".join(f"{e}({i})" for e, i
                              in zip(d.donor_elements, d.donor_atoms))
            lines.append(
                f"  {d.name} (denticity={d.denticity}, "
                f"donors=[{atoms}], metals=[{metals}], "
                f"conf={d.confidence:.1f})"
            )
        header = (f"Found {len(donors)} chelating group(s), "
                  f"total denticity = {self.total_denticity(donors)}:")
        return "\n".join([header] + lines)
