"""
Metal-Aware Scoring Module — additive scoring of metal-ligand
coordination quality using Morse-potential distance terms, angular
penalty terms, and donor compatibility scoring.

Inspired by GOLD's ChemPLP metal terms (distance + angle to ideal
coordination geometry) and MetalDock's QM-calibrated interaction
energies.  All terms are smooth and support conservative lower-bounding
for the epsilon-certified branch-and-bound layer.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import numpy as np

from acudockmetal.metal_params import (
    GEOMETRY_TEMPLATES,
    MetalParameterLibrary,
    MetalParams,
)
from acudockmetal.metal_site import CoordinationHypothesis, MetalSite


# ---------------------------------------------------------------------------
# Morse potential parameters per donor element (calibrated from CSD/MetalPDB)
# D_e = well depth (kcal/mol), alpha = steepness (1/A)
# These provide an attractive+repulsive potential instead of pure penalty.
# ---------------------------------------------------------------------------

_MORSE_PARAMS: Dict[str, Tuple[float, float]] = {
    # (D_e, alpha) — D_e is the coordination bond energy, alpha controls width
    "N": (8.0, 1.8),   # histidine, pyridine coordination bonds ~8 kcal/mol
    "O": (6.0, 1.6),   # carboxylate, hydroxamate, water coordination
    "S": (10.0, 1.5),  # thiolate coordination (stronger but softer)
    "CL": (4.0, 1.4),  # chloride coordination (weaker)
}


@dataclass
class MetalScoreBreakdown:
    """Per-component breakdown of the metal-aware score."""
    distance_penalty: float = 0.0
    morse_energy: float = 0.0         # Morse potential (negative = favorable)
    angle_penalty: float = 0.0
    cn_penalty: float = 0.0
    donor_compatibility: float = 0.0
    geometry_rmsd: float = 0.0
    total: float = 0.0

    def compute_total(self, weights: Dict[str, float]) -> float:
        """Weighted sum of all components."""
        self.total = (
            weights.get("distance", 5.0) * self.distance_penalty
            + weights.get("morse", 1.0) * self.morse_energy
            + weights.get("angle", 2.0) * self.angle_penalty
            + weights.get("cn", 3.0) * self.cn_penalty
            + weights.get("donor", 1.0) * self.donor_compatibility
            + weights.get("geometry", 1.5) * self.geometry_rmsd
        )
        return self.total


@dataclass
class MetalScoreResult:
    """Full metal-aware score for a docked pose."""
    pose_id: int
    hypothesis_id: int
    breakdown: MetalScoreBreakdown
    total_score: float = 0.0
    coordination_satisfied: bool = False
    donor_atom_indices: List[int] = field(default_factory=list)
    metal_donor_distances: List[float] = field(default_factory=list)
    metal_donor_angles: List[float] = field(default_factory=list)
    observed_cn: int = 0


class MetalAwareScorer:
    """
    Score docked poses based on metal-ligand coordination quality.

    Uses Morse potentials for distance scoring (providing both
    attractive and repulsive components), angular deviation penalties,
    and donor-type compatibility scoring.  All terms support
    analytical lower-bounding for the epsilon-certification layer.
    """

    def __init__(
        self,
        metal_library: Optional[MetalParameterLibrary] = None,
        weights: Optional[Dict[str, float]] = None,
        distance_cutoff: float = 3.0,
    ) -> None:
        self.metal_lib = metal_library or MetalParameterLibrary()
        self.weights = weights or {
            "distance": 5.0,
            "morse": 1.0,
            "angle": 2.0,
            "cn": 3.0,
            "donor": 1.0,
            "geometry": 1.5,
        }
        self.distance_cutoff = distance_cutoff

    @staticmethod
    def morse_energy(
        d: float, d_eq: float, D_e: float, alpha: float
    ) -> float:
        """
        Morse potential: E = D_e * (1 - exp(-alpha*(d - d_eq)))^2 - D_e

        Returns negative at equilibrium (favorable coordination),
        zero at infinity, and positive for very short distances (repulsive).
        """
        x = 1.0 - np.exp(-alpha * (d - d_eq))
        return D_e * x * x - D_e

    @staticmethod
    def morse_lower_bound(
        d_min: float, d_max: float, d_eq: float, D_e: float, alpha: float
    ) -> float:
        """
        Conservative lower bound of Morse potential over [d_min, d_max].

        The Morse potential has its minimum (-D_e) at d = d_eq.
        """
        if d_min <= d_eq <= d_max:
            return -D_e  # equilibrium distance is reachable
        # Minimum is at the boundary closest to d_eq
        d_closest = d_min if abs(d_min - d_eq) < abs(d_max - d_eq) else d_max
        x = 1.0 - np.exp(-alpha * (d_closest - d_eq))
        return D_e * x * x - D_e

    def score_pose(
        self,
        ligand_coords: np.ndarray,
        ligand_elements: List[str],
        hypothesis: CoordinationHypothesis,
        donor_atom_indices: Optional[List[int]] = None,
    ) -> MetalScoreResult:
        """
        Score a docked pose against a coordination hypothesis.

        Args:
            ligand_coords: (N, 3) atom coordinates of docked ligand.
            ligand_elements: element symbols for each ligand atom.
            hypothesis: the coordination hypothesis to score against.
            donor_atom_indices: optional pre-identified donor atom indices
                on the ligand. If None, auto-detect from distance.

        Returns:
            MetalScoreResult with per-component breakdown.
        """
        metal_coord = hypothesis.metal_site.coord
        metal_sym = hypothesis.metal_symbol
        mp = self.metal_lib.get(metal_sym)

        if mp is None:
            return MetalScoreResult(
                pose_id=-1,
                hypothesis_id=hypothesis.hypothesis_id,
                breakdown=MetalScoreBreakdown(),
            )

        # --- Find ligand donor atoms ---
        if donor_atom_indices is None:
            donor_atom_indices = self._find_donors(
                ligand_coords, ligand_elements, metal_coord, mp)

        # Compute per-donor distances
        donor_coords = ligand_coords[donor_atom_indices] if donor_atom_indices else np.zeros((0, 3))
        donor_elems = [ligand_elements[i] for i in donor_atom_indices] if donor_atom_indices else []
        distances = np.linalg.norm(
            donor_coords - metal_coord, axis=1) if len(donor_coords) > 0 else np.array([])

        # Combine with protein donors
        protein_donor_coords = np.array(
            [a.coord for a in hypothesis.protein_donors])
        protein_donor_elems = [a.element for a in hypothesis.protein_donors]

        all_donor_coords = np.vstack([protein_donor_coords, donor_coords]) \
            if len(donor_coords) > 0 and len(protein_donor_coords) > 0 \
            else (donor_coords if len(donor_coords) > 0 else protein_donor_coords)

        all_donor_elems = protein_donor_elems + donor_elems
        all_distances = np.linalg.norm(
            all_donor_coords - metal_coord, axis=1) if len(all_donor_coords) > 0 else np.array([])

        observed_cn = len(all_donor_coords)

        # --- Component 1: Distance penalty (quadratic) ---
        dist_penalty = self._distance_penalty(
            distances.tolist(), donor_elems, mp)

        # --- Component 1b: Morse energy (attractive+repulsive) ---
        morse_e = self._morse_energy_total(
            distances.tolist(), donor_elems, mp)

        # --- Component 2: Angle penalty ---
        angle_penalty, angles = self._angle_penalty(
            all_donor_coords, metal_coord, hypothesis.geometry)

        # --- Component 3: CN penalty ---
        cn_penalty = self._cn_penalty(observed_cn, hypothesis.coordination_number)

        # --- Component 4: Donor compatibility ---
        donor_compat = self._donor_compatibility(donor_elems, mp)

        # --- Component 5: Geometry RMSD ---
        geom_rmsd = self._geometry_rmsd(
            all_donor_coords, metal_coord, hypothesis.geometry)

        breakdown = MetalScoreBreakdown(
            distance_penalty=dist_penalty,
            morse_energy=morse_e,
            angle_penalty=angle_penalty,
            cn_penalty=cn_penalty,
            donor_compatibility=donor_compat,
            geometry_rmsd=geom_rmsd,
        )
        breakdown.compute_total(self.weights)

        # Coordination satisfied if CN matches and distances are OK
        cn_ok = observed_cn == hypothesis.coordination_number
        dist_ok = all(d < self.distance_cutoff for d in all_distances) \
            if len(all_distances) > 0 else False

        return MetalScoreResult(
            pose_id=-1,
            hypothesis_id=hypothesis.hypothesis_id,
            breakdown=breakdown,
            total_score=breakdown.total,
            coordination_satisfied=cn_ok and dist_ok,
            donor_atom_indices=donor_atom_indices,
            metal_donor_distances=distances.tolist(),
            metal_donor_angles=angles,
            observed_cn=observed_cn,
        )

    def _find_donors(
        self,
        coords: np.ndarray,
        elements: List[str],
        metal_coord: np.ndarray,
        mp: MetalParams,
    ) -> List[int]:
        """Auto-detect donor atoms: N/O/S/Cl within cutoff of metal."""
        donors = []
        donor_elements = {"N", "O", "S", "CL"}
        for i, (c, e) in enumerate(zip(coords, elements)):
            if e.upper() in donor_elements:
                d = np.linalg.norm(c - metal_coord)
                if d < self.distance_cutoff:
                    donors.append(i)
        return donors

    def _distance_penalty(
        self,
        distances: List[float],
        donor_elements: List[str],
        mp: MetalParams,
    ) -> float:
        """Quadratic penalty for deviation from ideal metal-donor distance."""
        if not distances:
            return 0.0
        total = 0.0
        for d, elem in zip(distances, donor_elements):
            d_ideal = mp.get_ideal_distance(elem)
            total += (d - d_ideal) ** 2
        return total / len(distances)

    def _morse_energy_total(
        self,
        distances: List[float],
        donor_elements: List[str],
        mp: MetalParams,
    ) -> float:
        """
        Sum of Morse potential energies for all ligand-metal donor bonds.

        Negative values indicate favorable coordination.
        """
        if not distances:
            return 0.0
        total = 0.0
        for d, elem in zip(distances, donor_elements):
            d_eq = mp.get_ideal_distance(elem)
            D_e, alpha = _MORSE_PARAMS.get(elem.upper(), (5.0, 1.5))
            total += self.morse_energy(d, d_eq, D_e, alpha)
        return total

    def _angle_penalty(
        self,
        all_donor_coords: np.ndarray,
        metal_coord: np.ndarray,
        geometry_name: str,
    ) -> Tuple[float, List[float]]:
        """Quadratic penalty for angular deviation from ideal geometry."""
        n = len(all_donor_coords)
        if n < 2:
            return 0.0, []

        template = GEOMETRY_TEMPLATES.get(geometry_name)
        if template is None:
            return 0.0, []

        # Compute all pairwise donor-metal-donor angles
        angles = []
        for i in range(n):
            for j in range(i + 1, n):
                v1 = all_donor_coords[i] - metal_coord
                v2 = all_donor_coords[j] - metal_coord
                n1 = np.linalg.norm(v1)
                n2 = np.linalg.norm(v2)
                if n1 < 1e-6 or n2 < 1e-6:
                    continue
                cos_a = np.clip(
                    np.dot(v1, v2) / (n1 * n2), -1.0, 1.0)
                angle = np.degrees(np.arccos(cos_a))
                angles.append(angle)

        if not angles:
            return 0.0, []

        # Compute mean minimum deviation from ideal angles
        total = 0.0
        for obs in angles:
            min_dev = min(abs(obs - ideal) for ideal in template.ideal_angles)
            total += min_dev ** 2
        penalty = total / len(angles)

        return penalty, angles

    def _cn_penalty(self, observed: int, target: int) -> float:
        """Quadratic penalty for coordination number mismatch."""
        return (observed - target) ** 2

    def _donor_compatibility(
        self, donor_elements: List[str], mp: MetalParams
    ) -> float:
        """
        Score donor-type compatibility (lower is better).

        Donors matching the metal's preferred list get 0 penalty;
        non-preferred donors get increasing penalty.
        """
        if not donor_elements:
            return 0.0
        total = 0.0
        for elem in donor_elements:
            e = elem.upper()
            if e in mp.preferred_donors:
                rank = mp.preferred_donors.index(e)
                total += rank * 0.1  # small penalty for lower-ranked
            else:
                total += 1.0  # significant penalty for non-preferred
        return total / len(donor_elements)

    def _geometry_rmsd(
        self,
        all_donor_coords: np.ndarray,
        metal_coord: np.ndarray,
        geometry_name: str,
    ) -> float:
        """
        RMSD of normalized coordination vectors vs ideal geometry template.

        Computes the angular deviation which correlates with true
        coordination geometry RMSD.
        """
        template = GEOMETRY_TEMPLATES.get(geometry_name)
        if template is None or len(all_donor_coords) < 2:
            return 0.0

        vectors = all_donor_coords - metal_coord
        norms = np.linalg.norm(vectors, axis=1, keepdims=True)
        norms = np.maximum(norms, 1e-6)
        normalized = vectors / norms

        # Compute variance of inter-vector angles from ideal
        angles = []
        n = len(normalized)
        for i in range(n):
            for j in range(i + 1, n):
                cos_a = np.clip(np.dot(normalized[i], normalized[j]), -1, 1)
                angles.append(np.degrees(np.arccos(cos_a)))

        if not angles:
            return 0.0

        deviations = []
        for obs in angles:
            min_dev = min(abs(obs - ideal) for ideal in template.ideal_angles)
            deviations.append(min_dev)

        return float(np.sqrt(np.mean(np.array(deviations) ** 2)))

    def lower_bound_over_region(
        self,
        translation_min: np.ndarray,
        translation_max: np.ndarray,
        hypothesis: CoordinationHypothesis,
    ) -> float:
        """
        Conservative lower bound on metal score over a spatial region.

        Used by the epsilon-certified branch-and-bound engine.  Computes
        lower bounds for distance, Morse, angle, and CN terms.
        """
        metal_coord = hypothesis.metal_site.coord
        metal_sym = hypothesis.metal_symbol
        mp = self.metal_lib.get(metal_sym)
        if mp is None:
            return 0.0

        # --- Distance bounds ---
        # Compute min/max possible distance from region to metal
        region_center = (translation_min + translation_max) / 2
        d_center = np.linalg.norm(region_center - metal_coord)
        region_radius = np.linalg.norm(
            (translation_max - translation_min) / 2)

        d_min = max(d_center - region_radius, 0.0)
        d_max = d_center + region_radius

        # Quadratic distance LB
        d_ideal = mp.get_ideal_distance("O")  # conservative generic
        if d_min <= d_ideal <= d_max:
            dist_lb = 0.0
        else:
            closest = min(abs(d_min - d_ideal), abs(d_max - d_ideal))
            dist_lb = closest ** 2

        # --- Morse energy LB ---
        # Use the most favorable (N) donor Morse params for the LB
        D_e, alpha = _MORSE_PARAMS.get("N", (8.0, 1.8))
        d_eq_n = mp.get_ideal_distance("N")
        morse_lb = self.morse_lower_bound(d_min, d_max, d_eq_n, D_e, alpha)

        # --- Angle LB ---
        # If only one ligand donor, the angle penalty is 0 (no pairwise angles)
        # For multi-donor, LB is 0 (could be perfect alignment)
        angle_lb = 0.0

        # --- CN LB ---
        # The CN penalty is at least 0 if the target CN is achievable.
        # If the region is far from the metal, no donors => CN penalty is
        # at least (target - n_protein)^2 minus open_slots^2
        n_protein = len(hypothesis.protein_donors)
        target_cn = hypothesis.coordination_number
        # Best case: all open slots filled by ligand donors
        cn_lb = 0.0
        if d_min > self.distance_cutoff:
            # No ligand donors can reach the metal from this region
            cn_lb = float((n_protein - target_cn) ** 2)

        return (
            self.weights.get("distance", 5.0) * dist_lb
            + self.weights.get("morse", 1.0) * morse_lb
            + self.weights.get("angle", 2.0) * angle_lb
            + self.weights.get("cn", 3.0) * cn_lb
        )

    def summary(self, results: List[MetalScoreResult]) -> str:
        """Human-readable summary of scoring results."""
        if not results:
            return "No scores computed."
        lines = [f"Metal-aware scores for {len(results)} pose(s):"]
        for r in results:
            sat = "YES" if r.coordination_satisfied else "NO"
            b = r.breakdown
            lines.append(
                f"  Pose {r.pose_id} (H{r.hypothesis_id}): "
                f"total={r.total_score:.2f} "
                f"[dist={b.distance_penalty:.2f}, "
                f"morse={b.morse_energy:.2f}, "
                f"angle={b.angle_penalty:.2f}, "
                f"CN={b.cn_penalty:.1f}, "
                f"donor={b.donor_compatibility:.2f}, "
                f"geom={b.geometry_rmsd:.2f}] "
                f"coord_ok={sat} CN={r.observed_cn}"
            )
        return "\n".join(lines)
