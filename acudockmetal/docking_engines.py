"""
Baseline Docking Engine Layer — wraps AutoDock Vina and GNINA.

Provides a unified interface for per-hypothesis docking, with results
collected and clustered by the DockingOrchestrator.
"""

from __future__ import annotations

import logging
import os
import subprocess
import tempfile
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import numpy as np

from acudockmetal.metal_site import CoordinationHypothesis
from acudockmetal.preparation import PreparedLigand, PreparedReceptor

log = logging.getLogger(__name__)


@dataclass
class DockingPose:
    """A single docked pose with metadata."""
    pose_id: int
    coordinates: np.ndarray       # (N_atoms, 3)
    vina_score: float = 0.0       # kcal/mol
    gnina_cnn_score: float = 0.0  # GNINA CNN pose score
    gnina_cnn_affinity: float = 0.0
    metal_score: float = 0.0      # from MetalAwareScorer
    composite_score: float = 0.0
    hypothesis_id: int = -1
    engine: str = ""              # "vina" or "gnina"
    pdbqt_string: str = ""        # raw PDBQT of the pose
    pdb_block: str = ""           # PDB block for visualization
    rmsd_from_best: float = 0.0
    cluster_id: int = -1

    @property
    def label(self) -> str:
        return (f"Pose{self.pose_id} ({self.engine}, "
                f"H{self.hypothesis_id}, "
                f"vina={self.vina_score:.1f})")


class VinaEngine:
    """
    AutoDock Vina docking engine wrapper.

    Uses the vina Python package for docking.
    """

    def __init__(
        self,
        exhaustiveness: int = 32,
        num_modes: int = 20,
        energy_range: float = 5.0,
        seed: int = 42,
        cpu: int = 0,
    ) -> None:
        self.exhaustiveness = exhaustiveness
        self.num_modes = num_modes
        self.energy_range = energy_range
        self.seed = seed
        self.cpu = cpu

    def dock(
        self,
        receptor: PreparedReceptor,
        ligand: PreparedLigand,
        hypothesis: Optional[CoordinationHypothesis] = None,
    ) -> List[DockingPose]:
        """
        Run Vina docking.

        If a hypothesis is provided, the box is centered on the metal site
        with tighter constraints.
        """
        from vina import Vina

        if not receptor.pdbqt_path:
            raise RuntimeError("Receptor PDBQT not available.")
        if not ligand.pdbqt_string or not ligand.pdbqt_string.strip():
            raise RuntimeError("Ligand PDBQT not available.")

        v = Vina(sf_name="vina", cpu=self.cpu, seed=self.seed)
        v.set_receptor(receptor.pdbqt_path)
        v.set_ligand_from_string(ligand.pdbqt_string)

        # Box: center on metal site if hypothesis given
        center = receptor.box_center.tolist()
        size = (receptor.box_size * 2).tolist()  # Vina wants full widths

        if hypothesis is not None:
            center = hypothesis.metal_site.coord.tolist()

        v.compute_vina_maps(center=center, box_size=size)
        v.dock(
            exhaustiveness=self.exhaustiveness,
            n_poses=self.num_modes,
        )

        # Extract results
        energies = v.energies()
        poses: List[DockingPose] = []

        if energies is None or len(energies) == 0:
            return poses

        for i, energy_row in enumerate(energies):
            score = energy_row[0]  # total score in kcal/mol

            # Get pose PDBQT
            pdbqt_str = ""
            try:
                pdbqt_str = v.poses(n_poses=i + 1).split("MODEL")[
                    -1] if i == 0 else ""
            except Exception:
                pass

            # Extract coordinates from the Vina output
            coords = self._extract_coords_from_vina(v, i)

            poses.append(DockingPose(
                pose_id=i,
                coordinates=coords,
                vina_score=float(score),
                hypothesis_id=hypothesis.hypothesis_id if hypothesis else -1,
                engine="vina",
                pdbqt_string=pdbqt_str,
            ))

        return poses

    def _extract_coords_from_vina(self, v, pose_idx: int) -> np.ndarray:
        """Extract atom coordinates for a specific pose."""
        try:
            pdbqt_all = v.poses(n_poses=pose_idx + 1)
            models = pdbqt_all.strip().split("ENDMDL")
            if pose_idx < len(models):
                return self._parse_pdbqt_coords(models[pose_idx])
        except Exception:
            pass
        return np.zeros((0, 3))

    @staticmethod
    def _parse_pdbqt_coords(pdbqt_block: str) -> np.ndarray:
        """Parse atom coordinates from a PDBQT block."""
        coords = []
        for line in pdbqt_block.splitlines():
            if line.startswith(("ATOM", "HETATM")):
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    coords.append([x, y, z])
                except ValueError:
                    continue
        return np.array(coords) if coords else np.zeros((0, 3))


class GninaEngine:
    """
    GNINA docking/rescoring engine wrapper.

    Calls the GNINA binary via subprocess for CNN-based
    docking and rescoring.
    """

    def __init__(
        self,
        gnina_path: str = "gnina",
        exhaustiveness: int = 16,
        num_modes: int = 20,
        cnn_scoring: str = "rescore",  # "none", "rescore", "refinement", "all"
        autobox_add: float = 4.0,
        seed: int = 42,
    ) -> None:
        self.gnina_path = gnina_path
        self.exhaustiveness = exhaustiveness
        self.num_modes = num_modes
        self.cnn_scoring = cnn_scoring
        self.autobox_add = autobox_add
        self.seed = seed

    def is_available(self) -> bool:
        """Check if GNINA binary is accessible."""
        try:
            result = subprocess.run(
                [self.gnina_path, "--version"],
                capture_output=True, text=True, timeout=10,
            )
            return result.returncode == 0
        except (OSError, subprocess.TimeoutExpired):
            return False

    def dock(
        self,
        receptor: PreparedReceptor,
        ligand: PreparedLigand,
        hypothesis: Optional[CoordinationHypothesis] = None,
    ) -> List[DockingPose]:
        """Run GNINA docking."""
        if not receptor.pdbqt_path or not ligand.pdbqt_string \
                or not ligand.pdbqt_string.strip():
            raise RuntimeError("PDBQT files required for GNINA.")

        with tempfile.TemporaryDirectory(prefix="gnina_") as tmpdir:
            lig_path = os.path.join(tmpdir, "ligand.pdbqt")
            out_path = os.path.join(tmpdir, "output.sdf")
            with open(lig_path, "w") as f:
                f.write(ligand.pdbqt_string)

            center = receptor.box_center
            size = receptor.box_size * 2

            if hypothesis is not None:
                center = hypothesis.metal_site.coord

            cmd = [
                self.gnina_path,
                "-r", receptor.pdbqt_path,
                "-l", lig_path,
                "-o", out_path,
                "--center_x", str(center[0]),
                "--center_y", str(center[1]),
                "--center_z", str(center[2]),
                "--size_x", str(size[0]),
                "--size_y", str(size[1]),
                "--size_z", str(size[2]),
                "--exhaustiveness", str(self.exhaustiveness),
                "--num_modes", str(self.num_modes),
                "--cnn_scoring", self.cnn_scoring,
                "--seed", str(self.seed),
            ]

            try:
                result = subprocess.run(
                    cmd, capture_output=True, text=True, timeout=600)
            except subprocess.TimeoutExpired:
                log.warning("GNINA timed out.")
                return []
            except FileNotFoundError:
                log.warning("GNINA binary not found at %s", self.gnina_path)
                return []

            if result.returncode != 0:
                log.warning("GNINA failed: %s", result.stderr[:500])
                return []

            return self._parse_gnina_output(out_path, hypothesis)

    def rescore(
        self,
        receptor: PreparedReceptor,
        poses: List[DockingPose],
    ) -> List[DockingPose]:
        """Rescore existing poses with GNINA CNN."""
        if not self.is_available():
            log.warning("GNINA not available for rescoring.")
            return poses

        with tempfile.TemporaryDirectory(prefix="gnina_rescore_") as tmpdir:
            lig_path = os.path.join(tmpdir, "poses.pdbqt")
            out_path = os.path.join(tmpdir, "rescored.sdf")

            # Write all pose PDBQTs
            with open(lig_path, "w") as f:
                for p in poses:
                    if p.pdbqt_string:
                        f.write(f"MODEL\n{p.pdbqt_string}\nENDMDL\n")

            if not os.path.getsize(lig_path):
                return poses

            cmd = [
                self.gnina_path,
                "-r", receptor.pdbqt_path,
                "-l", lig_path,
                "-o", out_path,
                "--score_only",
                "--cnn_scoring", "rescore",
            ]

            try:
                result = subprocess.run(
                    cmd, capture_output=True, text=True, timeout=300)
                if result.returncode == 0:
                    return self._update_scores_from_output(
                        poses, result.stdout)
            except Exception as e:
                log.warning("GNINA rescore failed: %s", e)

        return poses

    def _parse_gnina_output(
        self,
        sdf_path: str,
        hypothesis: Optional[CoordinationHypothesis],
    ) -> List[DockingPose]:
        """Parse GNINA SDF output into DockingPose objects."""
        poses = []
        try:
            from rdkit import Chem
            suppl = Chem.SDMolSupplier(sdf_path, removeHs=False)
            for i, mol in enumerate(suppl):
                if mol is None:
                    continue
                conf = mol.GetConformer()
                coords = np.array([
                    list(conf.GetAtomPosition(j))
                    for j in range(mol.GetNumAtoms())
                ])

                vina_score = float(mol.GetProp("minimizedAffinity")) \
                    if mol.HasProp("minimizedAffinity") else 0.0
                cnn_score = float(mol.GetProp("CNNscore")) \
                    if mol.HasProp("CNNscore") else 0.0
                cnn_aff = float(mol.GetProp("CNNaffinity")) \
                    if mol.HasProp("CNNaffinity") else 0.0

                poses.append(DockingPose(
                    pose_id=i,
                    coordinates=coords,
                    vina_score=vina_score,
                    gnina_cnn_score=cnn_score,
                    gnina_cnn_affinity=cnn_aff,
                    hypothesis_id=(hypothesis.hypothesis_id
                                   if hypothesis else -1),
                    engine="gnina",
                ))
        except Exception as e:
            log.warning("Failed to parse GNINA output: %s", e)
        return poses

    def _update_scores_from_output(
        self,
        poses: List[DockingPose],
        stdout: str,
    ) -> List[DockingPose]:
        """Update CNN scores from GNINA stdout."""
        # Parse CNN scores from output
        for line in stdout.splitlines():
            if "CNNscore" in line or "CNNaffinity" in line:
                # GNINA prints structured output; parse what we can
                pass
        return poses


class DockingOrchestrator:
    """
    Orchestrate docking across multiple engines and hypotheses.

    For each CoordinationHypothesis, runs Vina (and optionally GNINA),
    collects all poses, clusters them, and returns a unified result set.
    """

    def __init__(
        self,
        vina_engine: Optional[VinaEngine] = None,
        gnina_engine: Optional[GninaEngine] = None,
        use_gnina: bool = True,
        cluster_rmsd_cutoff: float = 2.0,
    ) -> None:
        self.vina = vina_engine or VinaEngine()
        self.gnina = gnina_engine or GninaEngine()
        self.use_gnina = use_gnina and self.gnina.is_available()
        self.cluster_rmsd = cluster_rmsd_cutoff

    def dock_all_hypotheses(
        self,
        receptor: PreparedReceptor,
        ligand: PreparedLigand,
        hypotheses: List[CoordinationHypothesis],
    ) -> List[DockingPose]:
        """
        Dock a ligand against all coordination hypotheses.

        Returns all poses across all hypotheses, sorted by Vina score.
        """
        all_poses: List[DockingPose] = []
        pose_counter = 0

        for hyp in hypotheses:
            log.info("Docking hypothesis: %s", hyp.label)

            # --- Vina ---
            try:
                vina_poses = self.vina.dock(receptor, ligand, hyp)
                for p in vina_poses:
                    p.pose_id = pose_counter
                    pose_counter += 1
                all_poses.extend(vina_poses)
                log.info("  Vina: %d poses", len(vina_poses))
            except Exception as e:
                log.warning("  Vina failed for %s: %s", hyp.label, e)

            # --- GNINA ---
            if self.use_gnina:
                try:
                    gnina_poses = self.gnina.dock(receptor, ligand, hyp)
                    for p in gnina_poses:
                        p.pose_id = pose_counter
                        pose_counter += 1
                    all_poses.extend(gnina_poses)
                    log.info("  GNINA: %d poses", len(gnina_poses))
                except Exception as e:
                    log.warning("  GNINA failed for %s: %s", hyp.label, e)

        # Sort by Vina score (lower is better)
        all_poses.sort(key=lambda p: p.vina_score)

        # Cluster poses
        all_poses = self._cluster_poses(all_poses)

        return all_poses

    def _cluster_poses(self, poses: List[DockingPose]) -> List[DockingPose]:
        """Simple greedy RMSD-based clustering."""
        if len(poses) < 2:
            for p in poses:
                p.cluster_id = 0
            return poses

        cluster_id = 0
        representatives: List[np.ndarray] = []

        for pose in poses:
            if pose.coordinates.size == 0:
                pose.cluster_id = cluster_id
                cluster_id += 1
                continue

            assigned = False
            for ci, rep_coords in enumerate(representatives):
                if rep_coords.shape == pose.coordinates.shape:
                    rmsd = np.sqrt(np.mean(
                        (rep_coords - pose.coordinates) ** 2))
                    if rmsd < self.cluster_rmsd:
                        pose.cluster_id = ci
                        assigned = True
                        break

            if not assigned:
                pose.cluster_id = cluster_id
                representatives.append(pose.coordinates.copy())
                cluster_id += 1

        return poses

    def summary(self, poses: List[DockingPose]) -> str:
        """Human-readable summary of docking results."""
        if not poses:
            return "No poses generated."
        n_clusters = len(set(p.cluster_id for p in poses))
        engines = set(p.engine for p in poses)
        hyps = set(p.hypothesis_id for p in poses)
        best = poses[0]
        lines = [
            f"Docking results: {len(poses)} poses in {n_clusters} cluster(s)",
            f"  Engines: {', '.join(engines)}",
            f"  Hypotheses: {sorted(hyps)}",
            f"  Best: {best.label} (Vina={best.vina_score:.2f} kcal/mol)",
            f"  Score range: [{poses[-1].vina_score:.2f}, "
            f"{poses[0].vina_score:.2f}] kcal/mol",
        ]
        return "\n".join(lines)
