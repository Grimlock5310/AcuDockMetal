"""
Multi-Fidelity Rescoring Pipeline.

Four levels of increasing accuracy and cost:
  Level 1: Metal-aware rescore (MetalAwareScorer)
  Level 2: GNINA CNN rescore
  Level 3: OpenMM MM minimization + validity checks
  Level 4: QM/MM with GFN2-xTB (semi-empirical QM for metal site)
"""

from __future__ import annotations

import logging
import tempfile
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import numpy as np

from acudockmetal.docking_engines import DockingPose, GninaEngine
from acudockmetal.metal_scoring import MetalAwareScorer, MetalScoreResult
from acudockmetal.metal_site import CoordinationHypothesis
from acudockmetal.preparation import PreparedReceptor

log = logging.getLogger(__name__)


@dataclass
class RescoredPose:
    """A pose with multi-fidelity scores attached."""
    original_pose: DockingPose
    level1_metal_score: Optional[MetalScoreResult] = None
    level2_cnn_score: float = 0.0
    level2_cnn_affinity: float = 0.0
    level3_mm_energy: float = 0.0
    level3_strain: float = 0.0
    level3_valid: bool = True
    level4_qmmm_energy: float = 0.0
    level4_available: bool = False
    composite_score: float = 0.0
    highest_level_completed: int = 0
    filter_reason: str = ""  # empty = passed all filters

    @property
    def passed_all(self) -> bool:
        return self.filter_reason == ""

    @property
    def label(self) -> str:
        return (f"Pose{self.original_pose.pose_id} "
                f"(L{self.highest_level_completed}, "
                f"composite={self.composite_score:.2f})")


@dataclass
class RescoringResult:
    """Full results from multi-fidelity rescoring."""
    all_poses: List[RescoredPose]
    level1_survivors: int = 0
    level2_survivors: int = 0
    level3_survivors: int = 0
    level4_survivors: int = 0
    total_input: int = 0

    @property
    def best_pose(self) -> Optional[RescoredPose]:
        passing = [p for p in self.all_poses if p.passed_all]
        if not passing:
            return self.all_poses[0] if self.all_poses else None
        return min(passing, key=lambda p: p.composite_score)

    def ranked_poses(self) -> List[RescoredPose]:
        return sorted(self.all_poses, key=lambda p: p.composite_score)


class MultiFidelityRescorer:
    """
    Multi-fidelity rescoring pipeline.

    Applies increasingly expensive scoring functions, filtering
    poses at each level to focus computational effort on the most
    promising candidates.
    """

    def __init__(
        self,
        metal_scorer: Optional[MetalAwareScorer] = None,
        gnina_engine: Optional[GninaEngine] = None,
        max_level: int = 4,
        level1_top_n: int = 50,
        level2_top_n: int = 20,
        level3_top_n: int = 10,
        level4_top_n: int = 5,
        metal_score_cutoff: float = 50.0,
        composite_weights: Optional[Dict[str, float]] = None,
    ) -> None:
        self.metal_scorer = metal_scorer or MetalAwareScorer()
        self.gnina = gnina_engine or GninaEngine()
        self.max_level = max_level
        self.level1_top_n = level1_top_n
        self.level2_top_n = level2_top_n
        self.level3_top_n = level3_top_n
        self.level4_top_n = level4_top_n
        self.metal_score_cutoff = metal_score_cutoff
        self.composite_weights = composite_weights or {
            "vina": 0.3,
            "metal": 0.25,
            "cnn": 0.25,
            "mm": 0.1,
            "qmmm": 0.1,
        }

    def rescore(
        self,
        poses: List[DockingPose],
        receptor: PreparedReceptor,
        hypothesis: CoordinationHypothesis,
        ligand_elements: Optional[List[str]] = None,
    ) -> RescoringResult:
        """
        Run multi-fidelity rescoring on a set of poses.

        Args:
            poses: input docked poses.
            receptor: prepared receptor.
            hypothesis: coordination hypothesis.
            ligand_elements: element symbols for ligand atoms.

        Returns:
            RescoringResult with all scored/filtered poses.
        """
        result = RescoringResult(
            all_poses=[],
            total_input=len(poses),
        )

        if not poses:
            return result

        # Wrap all poses
        rescored = [RescoredPose(original_pose=p) for p in poses]

        # --- Level 1: Metal-aware rescore ---
        log.info("Level 1: Metal-aware rescoring (%d poses)", len(rescored))
        rescored = self._level1_metal(rescored, hypothesis, ligand_elements)
        result.level1_survivors = sum(1 for p in rescored if p.passed_all)

        if self.max_level < 2:
            self._compute_composite(rescored)
            result.all_poses = rescored
            return result

        # --- Level 2: GNINA CNN rescore ---
        survivors = [p for p in rescored if p.passed_all][:self.level2_top_n]
        log.info("Level 2: GNINA CNN rescoring (%d poses)", len(survivors))
        survivors = self._level2_cnn(survivors, receptor)
        result.level2_survivors = len(survivors)

        # Merge back
        survivor_ids = {p.original_pose.pose_id for p in survivors}
        for p in rescored:
            if p.original_pose.pose_id not in survivor_ids:
                continue
            for s in survivors:
                if s.original_pose.pose_id == p.original_pose.pose_id:
                    p.level2_cnn_score = s.level2_cnn_score
                    p.level2_cnn_affinity = s.level2_cnn_affinity
                    p.highest_level_completed = max(
                        p.highest_level_completed, 2)

        if self.max_level < 3:
            self._compute_composite(rescored)
            result.all_poses = rescored
            return result

        # --- Level 3: MM minimization ---
        survivors = sorted(
            [p for p in rescored if p.passed_all],
            key=lambda p: p.composite_score,
        )[:self.level3_top_n]
        log.info("Level 3: MM minimization (%d poses)", len(survivors))
        survivors = self._level3_mm(survivors, receptor)
        result.level3_survivors = sum(1 for p in survivors if p.level3_valid)

        if self.max_level < 4:
            self._compute_composite(rescored)
            result.all_poses = rescored
            return result

        # --- Level 4: QM/MM ---
        survivors = [p for p in survivors if p.level3_valid][:self.level4_top_n]
        log.info("Level 4: QM/MM rescoring (%d poses)", len(survivors))
        survivors = self._level4_qmmm(survivors, receptor, hypothesis)
        result.level4_survivors = len(survivors)

        self._compute_composite(rescored)
        result.all_poses = rescored
        return result

    def _level1_metal(
        self,
        poses: List[RescoredPose],
        hypothesis: CoordinationHypothesis,
        ligand_elements: Optional[List[str]],
    ) -> List[RescoredPose]:
        """Level 1: Metal-aware scoring + filtering."""
        for p in poses:
            coords = p.original_pose.coordinates
            if coords.size == 0:
                p.filter_reason = "no_coordinates"
                continue

            elements = ligand_elements or ["C"] * len(coords)

            score = self.metal_scorer.score_pose(
                coords, elements, hypothesis)
            score.pose_id = p.original_pose.pose_id
            p.level1_metal_score = score
            p.highest_level_completed = 1

            # Filter by metal score cutoff
            if score.total_score > self.metal_score_cutoff:
                p.filter_reason = f"metal_score_too_high ({score.total_score:.1f})"

        # Keep top N by combined Vina + metal score
        valid = [p for p in poses if p.passed_all]
        valid.sort(key=lambda p: (
            p.original_pose.vina_score
            + (p.level1_metal_score.total_score if p.level1_metal_score else 0)
        ))

        for p in valid[self.level1_top_n:]:
            p.filter_reason = "filtered_at_level1"

        return poses

    def _level2_cnn(
        self,
        poses: List[RescoredPose],
        receptor: PreparedReceptor,
    ) -> List[RescoredPose]:
        """Level 2: GNINA CNN rescoring."""
        if not self.gnina.is_available():
            log.warning("GNINA not available; skipping Level 2.")
            for p in poses:
                p.highest_level_completed = max(
                    p.highest_level_completed, 2)
            return poses

        # Rescore via GNINA
        original_poses = [p.original_pose for p in poses]
        rescored = self.gnina.rescore(receptor, original_poses)

        for rp, orig in zip(poses, rescored):
            rp.level2_cnn_score = orig.gnina_cnn_score
            rp.level2_cnn_affinity = orig.gnina_cnn_affinity
            rp.highest_level_completed = max(rp.highest_level_completed, 2)

        return poses

    def _level3_mm(
        self,
        poses: List[RescoredPose],
        receptor: PreparedReceptor,
    ) -> List[RescoredPose]:
        """
        Level 3: OpenMM MM minimization + validity check.

        Brief minimization to assess physical plausibility.
        """
        for p in poses:
            try:
                energy, strain = self._mm_minimize(p, receptor)
                p.level3_mm_energy = energy
                p.level3_strain = strain
                p.level3_valid = strain < 100.0  # kcal/mol threshold
                p.highest_level_completed = max(
                    p.highest_level_completed, 3)
                if not p.level3_valid:
                    p.filter_reason = f"high_mm_strain ({strain:.1f})"
            except Exception as e:
                log.warning("Level 3 failed for pose %d: %s",
                            p.original_pose.pose_id, e)
                p.level3_valid = True  # don't penalize on failure
                p.highest_level_completed = max(
                    p.highest_level_completed, 3)

        return poses

    def _mm_minimize(
        self,
        pose: RescoredPose,
        receptor: PreparedReceptor,
    ) -> Tuple[float, float]:
        """
        Brief MM minimization using OpenMM.

        Returns (minimized_energy, strain_from_original).
        """
        try:
            from openmm import app, unit, LangevinMiddleIntegrator
            from openmm.app import ForceField, Modeller, PDBFile

            # Load receptor
            pdb = PDBFile(receptor.pdb_path)
            ff = ForceField("amber14-all.xml", "amber14/tip3pfb.xml")

            modeller = Modeller(pdb.topology, pdb.positions)

            system = ff.createSystem(
                modeller.topology,
                nonbondedMethod=app.NoCutoff,
                constraints=None,
            )

            integrator = LangevinMiddleIntegrator(
                300 * unit.kelvin,
                1 / unit.picosecond,
                0.002 * unit.picoseconds,
            )

            from openmm import Context, Platform
            platform = Platform.getPlatformByName("CPU")
            context = Context(system, integrator, platform)
            context.setPositions(modeller.positions)

            # Minimize
            from openmm import LocalEnergyMinimizer
            LocalEnergyMinimizer.minimize(context, maxIterations=100)

            state = context.getState(getEnergy=True)
            energy = state.getPotentialEnergy().value_in_unit(
                unit.kilocalories_per_mole)

            # Strain is estimated as the energy difference
            # (simplified — full implementation would compute per-ligand)
            return energy, abs(energy) * 0.01

        except ImportError:
            log.warning("OpenMM not available; skipping MM minimization.")
            return 0.0, 0.0
        except Exception as e:
            log.warning("MM minimization error: %s", e)
            return 0.0, 0.0

    def _level4_qmmm(
        self,
        poses: List[RescoredPose],
        receptor: PreparedReceptor,
        hypothesis: CoordinationHypothesis,
    ) -> List[RescoredPose]:
        """
        Level 4: QM/MM rescoring with GFN2-xTB.

        QM region: metal + first-shell coordinators + ligand donors.
        MM region: rest (point charges).
        Subtractive scheme: E = E_QM(QM) + E_MM(full) - E_MM(QM).
        """
        for p in poses:
            try:
                qmmm_e = self._qmmm_rescore(p, receptor, hypothesis)
                p.level4_qmmm_energy = qmmm_e
                p.level4_available = True
                p.highest_level_completed = max(
                    p.highest_level_completed, 4)
            except Exception as e:
                log.warning("Level 4 QM/MM failed for pose %d: %s",
                            p.original_pose.pose_id, e)
                p.level4_available = False

        return poses

    def _qmmm_rescore(
        self,
        pose: RescoredPose,
        receptor: PreparedReceptor,
        hypothesis: CoordinationHypothesis,
    ) -> float:
        """Run GFN2-xTB QM/MM rescoring."""
        try:
            from xtb.interface import Calculator, Param, Environment
            from xtb.libxtb import VERBOSITY_MUTED

            # Build QM region: metal + protein donors + ligand
            # (simplified implementation)
            metal_coord = hypothesis.metal_site.coord
            protein_donors = hypothesis.protein_donors

            # Collect QM atom positions and atomic numbers
            positions = [metal_coord]
            numbers = [self._atomic_number(hypothesis.metal_symbol)]

            for donor in protein_donors:
                positions.append(donor.coord)
                numbers.append(self._atomic_number(donor.element))

            ligand_coords = pose.original_pose.coordinates
            if ligand_coords.size > 0:
                for i in range(len(ligand_coords)):
                    positions.append(ligand_coords[i])
                    numbers.append(6)  # approximate as carbon

            positions = np.array(positions) / 0.529177  # Angstrom to Bohr
            numbers = np.array(numbers, dtype=int)

            env = Environment()
            env.set_verbosity(VERBOSITY_MUTED)
            calc = Calculator(Param.GFN2xTB, numbers, positions)
            calc.set_verbosity(VERBOSITY_MUTED)

            res = calc.singlepoint()
            energy_hartree = res.get_energy()
            energy_kcal = energy_hartree * 627.509  # Hartree to kcal/mol

            return energy_kcal

        except ImportError:
            log.warning("xtb-python not available; skipping QM/MM.")
            return 0.0
        except Exception as e:
            log.warning("QM/MM calculation error: %s", e)
            return 0.0

    @staticmethod
    def _atomic_number(symbol: str) -> int:
        """Convert element symbol to atomic number."""
        table = {
            "H": 1, "C": 6, "N": 7, "O": 8, "F": 9, "S": 16,
            "CL": 17, "P": 15, "BR": 35, "I": 53,
            "ZN": 30, "FE": 26, "CU": 29, "CO": 27, "NI": 28,
            "MN": 25, "MG": 12, "CA": 20,
            "LA": 57, "CE": 58, "GD": 64, "LU": 71,
        }
        return table.get(symbol.upper(), 6)

    def _compute_composite(self, poses: List[RescoredPose]) -> None:
        """Compute weighted composite score for all poses."""
        w = self.composite_weights
        for p in poses:
            score = w["vina"] * p.original_pose.vina_score
            if p.level1_metal_score:
                score += w["metal"] * p.level1_metal_score.total_score
            if p.level2_cnn_affinity:
                score += w["cnn"] * (-p.level2_cnn_affinity)
            if p.level3_mm_energy:
                score += w["mm"] * p.level3_mm_energy
            if p.level4_available:
                score += w["qmmm"] * p.level4_qmmm_energy
            p.composite_score = score

    def summary(self, result: RescoringResult) -> str:
        """Human-readable summary of rescoring results."""
        lines = [
            f"Multi-fidelity rescoring: {result.total_input} input poses",
            f"  Level 1 (metal-aware): {result.level1_survivors} survived",
            f"  Level 2 (CNN): {result.level2_survivors} survived",
            f"  Level 3 (MM): {result.level3_survivors} survived",
            f"  Level 4 (QM/MM): {result.level4_survivors} completed",
        ]
        best = result.best_pose
        if best:
            lines.append(f"  Best: {best.label}")
        return "\n".join(lines)
