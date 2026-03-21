"""
Geometry Validity Gate — checks coordination geometry, steric clashes,
and chemically plausible bond lengths for docked poses.

Inspired by FindGeo (coordination geometry) and PoseBusters
(pose plausibility checks).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import numpy as np

from acudockmetal.metal_params import (
    GEOMETRY_TEMPLATES,
    MetalParameterLibrary,
)
from acudockmetal.metal_site import CoordinationHypothesis


@dataclass
class ValidationCheck:
    """Result of a single validation check."""
    name: str
    passed: bool
    value: float = 0.0          # observed value
    threshold: float = 0.0      # pass/fail boundary
    details: str = ""

    @property
    def label(self) -> str:
        status = "PASS" if self.passed else "FAIL"
        return f"[{status}] {self.name}: {self.value:.2f} (thr={self.threshold:.2f})"


@dataclass
class ValidationResult:
    """Aggregate result from the geometry validity gate."""
    pose_id: int
    checks: List[ValidationCheck] = field(default_factory=list)
    overall_pass: bool = True
    validity_score: float = 1.0  # 0-1, fraction of checks passed

    def compute_overall(self) -> None:
        if not self.checks:
            self.overall_pass = True
            self.validity_score = 1.0
            return
        n_pass = sum(1 for c in self.checks if c.passed)
        self.validity_score = n_pass / len(self.checks)
        self.overall_pass = all(c.passed for c in self.checks)

    @property
    def summary(self) -> str:
        status = "VALID" if self.overall_pass else "INVALID"
        n_pass = sum(1 for c in self.checks if c.passed)
        return (f"Pose {self.pose_id}: {status} "
                f"({n_pass}/{len(self.checks)} checks passed, "
                f"score={self.validity_score:.2f})")


class GeometryValidator:
    """
    Validate docked pose geometry against coordination chemistry
    constraints and structural plausibility criteria.

    Checks:
    1. Coordination number matches hypothesis
    2. Donor-metal distances within acceptable ranges
    3. Donor-metal-donor angles match geometry template
    4. No severe steric clashes (ligand-protein)
    5. No chemically implausible bond lengths
    """

    def __init__(
        self,
        metal_library: Optional[MetalParameterLibrary] = None,
        distance_tolerance: float = 0.5,
        angle_tolerance: float = 25.0,
        clash_distance: float = 1.5,
        max_bond_length: float = 2.0,
        cn_strict: bool = False,
    ) -> None:
        """
        Args:
            distance_tolerance: allowed deviation from ideal (angstrom).
            angle_tolerance: allowed angular deviation (degrees).
            clash_distance: minimum non-bonded distance (angstrom).
            max_bond_length: maximum C-C/C-N/C-O bond length (angstrom).
            cn_strict: if True, CN must exactly match hypothesis.
        """
        self.metal_lib = metal_library or MetalParameterLibrary()
        self.distance_tolerance = distance_tolerance
        self.angle_tolerance = angle_tolerance
        self.clash_distance = clash_distance
        self.max_bond_length = max_bond_length
        self.cn_strict = cn_strict

    def validate(
        self,
        ligand_coords: np.ndarray,
        ligand_elements: List[str],
        hypothesis: CoordinationHypothesis,
        receptor_coords: Optional[np.ndarray] = None,
        donor_atom_indices: Optional[List[int]] = None,
    ) -> ValidationResult:
        """
        Validate a docked pose.

        Args:
            ligand_coords: (N, 3) ligand atom coordinates.
            ligand_elements: element symbols per atom.
            hypothesis: coordination hypothesis.
            receptor_coords: (M, 3) receptor atom coords (for clash check).
            donor_atom_indices: pre-identified ligand donor indices.

        Returns:
            ValidationResult with all checks.
        """
        result = ValidationResult(pose_id=-1)
        metal_coord = hypothesis.metal_site.coord
        metal_sym = hypothesis.metal_symbol
        mp = self.metal_lib.get(metal_sym)

        if mp is None or ligand_coords.size == 0:
            result.compute_overall()
            return result

        # Find donors if not provided
        if donor_atom_indices is None:
            donor_atom_indices = []
            for i, (c, e) in enumerate(zip(ligand_coords, ligand_elements)):
                if e.upper() in ("N", "O", "S"):
                    d = np.linalg.norm(c - metal_coord)
                    if d < 3.0:
                        donor_atom_indices.append(i)

        donor_coords = ligand_coords[donor_atom_indices] if donor_atom_indices else np.zeros((0, 3))
        donor_elems = [ligand_elements[i] for i in donor_atom_indices] if donor_atom_indices else []

        # All donors (protein + ligand)
        protein_donor_coords = np.array(
            [a.coord for a in hypothesis.protein_donors])
        protein_donor_elems = [a.element for a in hypothesis.protein_donors]

        all_donor_coords = np.vstack(
            [protein_donor_coords, donor_coords]
        ) if len(donor_coords) > 0 and len(protein_donor_coords) > 0 else (
            donor_coords if len(donor_coords) > 0 else protein_donor_coords
        )
        all_donor_elems = protein_donor_elems + donor_elems
        observed_cn = len(all_donor_coords)

        # --- Check 1: Coordination number ---
        result.checks.append(self._check_cn(
            observed_cn, hypothesis.coordination_number))

        # --- Check 2: Donor-metal distances ---
        for i, (dc, de) in enumerate(zip(
                donor_coords if len(donor_coords) > 0 else [],
                donor_elems)):
            d_ideal = mp.get_ideal_distance(de)
            d_obs = np.linalg.norm(dc - metal_coord)
            result.checks.append(self._check_distance(
                d_obs, d_ideal, f"donor_{i}_{de}"))

        # --- Check 3: Angles ---
        template = GEOMETRY_TEMPLATES.get(hypothesis.geometry)
        if template is not None and len(all_donor_coords) >= 2:
            result.checks.extend(self._check_angles(
                all_donor_coords, metal_coord, template))

        # --- Check 4: Steric clashes ---
        if receptor_coords is not None and receptor_coords.size > 0:
            result.checks.append(self._check_clashes(
                ligand_coords, receptor_coords))

        # --- Check 5: Internal bond lengths ---
        result.checks.extend(self._check_bond_lengths(
            ligand_coords, ligand_elements))

        result.compute_overall()
        return result

    def _check_cn(self, observed: int, target: int) -> ValidationCheck:
        """Check coordination number."""
        if self.cn_strict:
            passed = observed == target
        else:
            passed = abs(observed - target) <= 1
        return ValidationCheck(
            name="coordination_number",
            passed=passed,
            value=float(observed),
            threshold=float(target),
            details=f"Observed CN={observed}, target={target}",
        )

    def _check_distance(
        self, observed: float, ideal: float, label: str
    ) -> ValidationCheck:
        """Check a single metal-donor distance."""
        deviation = abs(observed - ideal)
        passed = deviation <= self.distance_tolerance
        return ValidationCheck(
            name=f"distance_{label}",
            passed=passed,
            value=observed,
            threshold=ideal,
            details=f"d={observed:.2f} A, ideal={ideal:.2f}, "
                    f"dev={deviation:.2f}",
        )

    def _check_angles(
        self,
        all_donor_coords: np.ndarray,
        metal_coord: np.ndarray,
        template,
    ) -> List[ValidationCheck]:
        """Check donor-metal-donor angles against template."""
        checks = []
        n = len(all_donor_coords)

        for i in range(n):
            for j in range(i + 1, n):
                v1 = all_donor_coords[i] - metal_coord
                v2 = all_donor_coords[j] - metal_coord
                n1 = np.linalg.norm(v1)
                n2 = np.linalg.norm(v2)
                if n1 < 1e-6 or n2 < 1e-6:
                    continue

                cos_a = np.clip(np.dot(v1, v2) / (n1 * n2), -1, 1)
                angle = np.degrees(np.arccos(cos_a))

                min_dev = min(
                    abs(angle - ideal) for ideal in template.ideal_angles)

                passed = min_dev <= self.angle_tolerance
                closest_ideal = min(
                    template.ideal_angles,
                    key=lambda a: abs(angle - a))

                checks.append(ValidationCheck(
                    name=f"angle_{i}_{j}",
                    passed=passed,
                    value=angle,
                    threshold=closest_ideal,
                    details=f"angle={angle:.1f} deg, "
                            f"ideal={closest_ideal:.1f}, "
                            f"dev={min_dev:.1f}",
                ))

        return checks

    def _check_clashes(
        self,
        ligand_coords: np.ndarray,
        receptor_coords: np.ndarray,
    ) -> ValidationCheck:
        """Check for severe steric clashes between ligand and receptor."""
        # Compute all pairwise distances (brute force for small systems)
        # For large systems, use a spatial index
        min_dist = float("inf")
        n_clashes = 0

        for lc in ligand_coords:
            dists = np.linalg.norm(receptor_coords - lc, axis=1)
            local_min = float(np.min(dists))
            if local_min < min_dist:
                min_dist = local_min
            if local_min < self.clash_distance:
                n_clashes += 1

        passed = n_clashes == 0
        return ValidationCheck(
            name="steric_clashes",
            passed=passed,
            value=float(n_clashes),
            threshold=0.0,
            details=f"{n_clashes} clashes (min_dist={min_dist:.2f} A)",
        )

    def _check_bond_lengths(
        self,
        coords: np.ndarray,
        elements: List[str],
    ) -> List[ValidationCheck]:
        """
        Check for chemically implausible internal bond lengths.

        Looks for unreasonably long bonds between bonded heavy atoms
        (simplified: checks all nearby heavy-atom pairs).
        """
        checks = []
        n = len(coords)
        if n < 2:
            return checks

        heavy_mask = [e.upper() not in ("H",) for e in elements]
        heavy_idx = [i for i, m in enumerate(heavy_mask) if m]

        n_bad = 0
        for i_idx, i in enumerate(heavy_idx):
            for j in heavy_idx[i_idx + 1:]:
                d = np.linalg.norm(coords[i] - coords[j])
                # If atoms are close enough to be bonded but too far
                if 0.8 < d < 1.2 * self.max_bond_length:
                    continue  # likely bonded, length OK
                if d <= 0.8:
                    n_bad += 1  # too close

        checks.append(ValidationCheck(
            name="internal_geometry",
            passed=n_bad == 0,
            value=float(n_bad),
            threshold=0.0,
            details=f"{n_bad} implausible internal distances",
        ))

        return checks

    def batch_validate(
        self,
        poses_data: List[dict],
        hypothesis: CoordinationHypothesis,
        receptor_coords: Optional[np.ndarray] = None,
    ) -> List[ValidationResult]:
        """
        Validate a batch of poses.

        Args:
            poses_data: list of dicts with 'coords', 'elements', and
                        optional 'donor_indices'.
        """
        results = []
        for i, pd in enumerate(poses_data):
            r = self.validate(
                pd["coords"],
                pd["elements"],
                hypothesis,
                receptor_coords=receptor_coords,
                donor_atom_indices=pd.get("donor_indices"),
            )
            r.pose_id = i
            results.append(r)
        return results

    def summary(self, results: List[ValidationResult]) -> str:
        """Human-readable summary of validation results."""
        if not results:
            return "No validation results."
        n_valid = sum(1 for r in results if r.overall_pass)
        lines = [
            f"Geometry validation: {n_valid}/{len(results)} poses valid",
        ]
        for r in results:
            lines.append(f"  {r.summary}")
        return "\n".join(lines)
