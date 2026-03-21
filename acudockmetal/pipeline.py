"""
Pipeline Orchestrator — ties all AcuDockMetal components into a single
end-to-end docking workflow.
"""

from __future__ import annotations

import logging
import time
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional

import numpy as np

from acudockmetal.metal_params import MetalParameterLibrary
from acudockmetal.chelator_detect import ChelatingGroupDetector
from acudockmetal.metal_site import (
    CoordinationHypothesis,
    CoordinationHypothesisGenerator,
    MetalSite,
    MetalSiteDetector,
)
from acudockmetal.preparation import (
    LigandPreparator,
    PreparedLigand,
    PreparedReceptor,
    ReceptorPreparator,
)
from acudockmetal.docking_engines import (
    DockingOrchestrator,
    DockingPose,
    GninaEngine,
    VinaEngine,
)
from acudockmetal.metal_scoring import MetalAwareScorer, MetalScoreResult
from acudockmetal.certified_dock import CertifiedDockingEngine, CertificationResult
from acudockmetal.rescoring import MultiFidelityRescorer, RescoringResult
from acudockmetal.validation import GeometryValidator, ValidationResult
from acudockmetal.visualization import DockingVisualizer

log = logging.getLogger(__name__)


@dataclass
class AcuDockConfig:
    """Global configuration for the AcuDockMetal pipeline."""
    # Docking
    exhaustiveness: int = 32
    num_poses: int = 20
    box_padding: float = 10.0

    # Certification
    epsilon: float = 1.0
    certification_time_limit: float = 300.0
    certification_node_limit: int = 500_000
    run_certification: bool = True

    # Rescoring
    max_rescore_level: int = 4
    qm_mm_enabled: bool = True
    gnina_enabled: bool = True

    # Metal
    metals_enabled: List[str] = field(
        default_factory=lambda: ["ZN", "FE", "CU", "CO", "NI", "MN",
                                  "MG", "CA", "LN"])
    max_hypotheses_per_site: int = 12

    # Preparation
    ph: float = 7.4
    enumerate_tautomers: bool = True
    max_tautomers: int = 5
    num_conformers: int = 10

    # Visualization
    viz_width: int = 800
    viz_height: int = 500


@dataclass
class PipelineResult:
    """Complete results from a pipeline run."""
    # Inputs
    receptor_path: str = ""
    ligand_source: str = ""
    config: Optional[AcuDockConfig] = None

    # Preparation
    prepared_receptor: Optional[PreparedReceptor] = None
    prepared_ligands: List[PreparedLigand] = field(default_factory=list)
    metal_sites: List[MetalSite] = field(default_factory=list)
    hypotheses: List[CoordinationHypothesis] = field(default_factory=list)

    # Docking
    docking_poses: List[DockingPose] = field(default_factory=list)

    # Scoring
    metal_scores: List[MetalScoreResult] = field(default_factory=list)
    validations: List[ValidationResult] = field(default_factory=list)

    # Certification
    certification: Optional[CertificationResult] = None

    # Rescoring
    rescoring: Optional[RescoringResult] = None

    # Timing
    timing: Dict[str, float] = field(default_factory=dict)

    @property
    def best_pose(self) -> Optional[DockingPose]:
        if self.docking_poses:
            return self.docking_poses[0]
        return None

    @property
    def summary(self) -> str:
        lines = [
            "=" * 60,
            "AcuDockMetal Pipeline Results",
            "=" * 60,
            f"Receptor: {self.receptor_path}",
            f"Ligand: {self.ligand_source}",
            f"Metal sites: {len(self.metal_sites)}",
            f"Hypotheses: {len(self.hypotheses)}",
            f"Docking poses: {len(self.docking_poses)}",
        ]
        if self.best_pose:
            lines.append(
                f"Best pose: Vina={self.best_pose.vina_score:.2f} kcal/mol")
        if self.certification:
            lines.append(
                f"Certification: {self.certification.certification_status} "
                f"(gap={self.certification.actual_gap:.3f})")
        if self.rescoring:
            best_r = self.rescoring.best_pose
            if best_r:
                lines.append(
                    f"Best rescored: {best_r.label}")

        # Timing
        total = sum(self.timing.values())
        lines.append(f"\nTiming (total={total:.1f}s):")
        for step, t in self.timing.items():
            lines.append(f"  {step}: {t:.1f}s")

        lines.append("=" * 60)
        return "\n".join(lines)


class AcuDockMetalPipeline:
    """
    End-to-end metal-aware certified docking pipeline.

    Steps:
    1. Prepare receptor (PDBFixer + metal site detection)
    2. Prepare ligand (RDKit + chelator detection + tautomers)
    3. Generate coordination hypotheses
    4. Baseline docking (Vina + GNINA) per hypothesis
    5. Metal-aware scoring
    6. Geometry validation
    7. (Optional) epsilon-certified audit
    8. Multi-fidelity rescoring
    9. Visualization
    """

    def __init__(self, config: Optional[AcuDockConfig] = None) -> None:
        self.config = config or AcuDockConfig()
        self.metal_lib = MetalParameterLibrary()

        # Initialize components
        self.receptor_prep = ReceptorPreparator(
            ph=self.config.ph,
            box_padding=self.config.box_padding,
        )
        self.ligand_prep = LigandPreparator(
            num_conformers=self.config.num_conformers,
            enumerate_tautomers=self.config.enumerate_tautomers,
            max_tautomers=self.config.max_tautomers,
        )
        self.site_detector = MetalSiteDetector(
            metal_library=self.metal_lib)
        self.hyp_generator = CoordinationHypothesisGenerator(
            metal_library=self.metal_lib,
            max_hypotheses_per_site=self.config.max_hypotheses_per_site,
        )
        self.vina = VinaEngine(
            exhaustiveness=self.config.exhaustiveness,
            num_modes=self.config.num_poses,
        )
        self.gnina = GninaEngine(
            num_modes=self.config.num_poses,
        )
        self.orchestrator = DockingOrchestrator(
            vina_engine=self.vina,
            gnina_engine=self.gnina,
            use_gnina=self.config.gnina_enabled,
        )
        self.metal_scorer = MetalAwareScorer(
            metal_library=self.metal_lib)
        self.validator = GeometryValidator(
            metal_library=self.metal_lib)
        self.certified_engine = CertifiedDockingEngine(
            metal_scorer=self.metal_scorer,
            epsilon=self.config.epsilon,
            time_limit=self.config.certification_time_limit,
            node_limit=self.config.certification_node_limit,
        )
        self.rescorer = MultiFidelityRescorer(
            metal_scorer=self.metal_scorer,
            gnina_engine=self.gnina,
            max_level=self.config.max_rescore_level,
        )
        self.visualizer = DockingVisualizer(
            width=self.config.viz_width,
            height=self.config.viz_height,
        )

    def run(
        self,
        receptor_path: str,
        ligand_input: str,
        output_dir: Optional[str] = None,
    ) -> PipelineResult:
        """
        Run the full docking pipeline.

        Args:
            receptor_path: path to receptor PDB file.
            ligand_input: SMILES string or path to SDF/MOL file.
            output_dir: directory for output files.

        Returns:
            PipelineResult with all results.
        """
        result = PipelineResult(
            receptor_path=receptor_path,
            ligand_source=ligand_input,
            config=self.config,
        )

        # --- Step 1: Prepare receptor ---
        t0 = time.time()
        log.info("Step 1: Preparing receptor...")
        try:
            receptor = self.receptor_prep.prepare(
                receptor_path, output_dir=output_dir)
            result.prepared_receptor = receptor
            result.metal_sites = receptor.metal_sites
        except Exception as e:
            log.error("Receptor preparation failed: %s", e)
            result.timing["preparation"] = time.time() - t0
            return result
        result.timing["receptor_prep"] = time.time() - t0

        # --- Step 2: Prepare ligand ---
        t0 = time.time()
        log.info("Step 2: Preparing ligand...")
        try:
            if ligand_input.endswith((".sdf", ".mol")):
                ligands = self.ligand_prep.prepare_from_sdf(
                    ligand_input, output_dir=output_dir)
            else:
                ligands = self.ligand_prep.prepare_from_smiles(
                    ligand_input, output_dir=output_dir)
            result.prepared_ligands = ligands
        except Exception as e:
            log.error("Ligand preparation failed: %s", e)
            result.timing["ligand_prep"] = time.time() - t0
            return result
        result.timing["ligand_prep"] = time.time() - t0

        if not ligands:
            log.error("No ligand variants prepared.")
            return result

        # --- Step 3: Generate hypotheses ---
        t0 = time.time()
        log.info("Step 3: Generating coordination hypotheses...")
        all_hypotheses: List[CoordinationHypothesis] = []
        for site in result.metal_sites:
            hyps = self.hyp_generator.generate(site)
            all_hypotheses.extend(hyps)
        result.hypotheses = all_hypotheses
        log.info("  Generated %d hypotheses for %d metal site(s)",
                 len(all_hypotheses), len(result.metal_sites))
        result.timing["hypothesis_gen"] = time.time() - t0

        if not all_hypotheses:
            log.warning("No hypotheses generated; running generic docking.")
            # Could still dock without metal-specific constraints

        # --- Step 4: Baseline docking ---
        t0 = time.time()
        log.info("Step 4: Baseline docking...")
        ligand = ligands[0]  # primary variant

        if all_hypotheses:
            poses = self.orchestrator.dock_all_hypotheses(
                receptor, ligand, all_hypotheses)
        else:
            # Fall back to standard docking
            try:
                poses = self.vina.dock(receptor, ligand)
            except Exception as e:
                log.error("Docking failed: %s", e)
                poses = []

        result.docking_poses = poses
        log.info("  Generated %d poses", len(poses))
        result.timing["docking"] = time.time() - t0

        if not poses:
            return result

        # --- Step 5: Metal-aware scoring ---
        t0 = time.time()
        log.info("Step 5: Metal-aware scoring...")
        if all_hypotheses:
            for pose in poses:
                hyp = self._get_hypothesis(
                    all_hypotheses, pose.hypothesis_id)
                if hyp is None:
                    continue
                # Get ligand elements
                elements = self._get_elements(ligand)
                score = self.metal_scorer.score_pose(
                    pose.coordinates, elements, hyp)
                score.pose_id = pose.pose_id
                result.metal_scores.append(score)
                pose.metal_score = score.total_score
                pose.composite_score = (
                    pose.vina_score + 0.3 * score.total_score)
        result.timing["metal_scoring"] = time.time() - t0

        # --- Step 6: Geometry validation ---
        t0 = time.time()
        log.info("Step 6: Geometry validation...")
        if all_hypotheses:
            for pose in poses[:50]:  # validate top 50
                hyp = self._get_hypothesis(
                    all_hypotheses, pose.hypothesis_id)
                if hyp is None:
                    continue
                elements = self._get_elements(ligand)
                val = self.validator.validate(
                    pose.coordinates, elements, hyp)
                val.pose_id = pose.pose_id
                result.validations.append(val)
        result.timing["validation"] = time.time() - t0

        # --- Step 7: ε-Certified audit ---
        if self.config.run_certification and all_hypotheses:
            t0 = time.time()
            log.info("Step 7: epsilon-certified docking audit...")
            best_hyp = all_hypotheses[0]
            cert = self.certified_engine.certified_dock(
                hypothesis=best_hyp,
                warm_start_poses=poses[:5],
                box_center=receptor.box_center,
                box_half_width=self.config.box_padding,
            )
            result.certification = cert
            log.info("  %s", cert.certification_status)
            result.timing["certification"] = time.time() - t0

        # --- Step 8: Multi-fidelity rescoring ---
        t0 = time.time()
        log.info("Step 8: Multi-fidelity rescoring...")
        if all_hypotheses:
            hyp = all_hypotheses[0]
            elements = self._get_elements(ligand)
            rescore = self.rescorer.rescore(
                poses, receptor, hyp, ligand_elements=elements)
            result.rescoring = rescore
        result.timing["rescoring"] = time.time() - t0

        # Re-sort poses by composite score
        result.docking_poses.sort(
            key=lambda p: p.composite_score)

        log.info("Pipeline complete.")
        return result

    def _get_hypothesis(
        self,
        hypotheses: List[CoordinationHypothesis],
        hyp_id: int,
    ) -> Optional[CoordinationHypothesis]:
        """Look up a hypothesis by ID."""
        for h in hypotheses:
            if h.hypothesis_id == hyp_id:
                return h
        return hypotheses[0] if hypotheses else None

    def _get_elements(self, ligand: PreparedLigand) -> List[str]:
        """Extract element symbols from a prepared ligand."""
        try:
            mol = ligand.mol
            return [atom.GetSymbol() for atom in mol.GetAtoms()]
        except Exception:
            return []

    # ------------------------------------------------------------------
    # Convenience visualization methods
    # ------------------------------------------------------------------

    def show_metal_site(self, result: PipelineResult, site_idx: int = 0):
        """Show 3D view of a metal site."""
        if not result.metal_sites or not result.prepared_receptor:
            log.warning("No metal sites to display.")
            return None
        site = result.metal_sites[site_idx]
        return self.visualizer.view_metal_site(
            result.prepared_receptor.pdb_path, site)

    def show_best_pose(self, result: PipelineResult):
        """Show 3D view of the best docked pose."""
        if not result.docking_poses or not result.prepared_receptor:
            log.warning("No poses to display.")
            return None
        hyp = result.hypotheses[0] if result.hypotheses else None
        return self.visualizer.view_docked_pose(
            result.prepared_receptor.pdb_path,
            result.docking_poses[0],
            hypothesis=hyp,
        )

    def plot_results(self, result: PipelineResult) -> Dict[str, Any]:
        """Generate all 2D analysis plots."""
        plots = {}
        if result.docking_poses:
            plots["score_dist"] = self.visualizer.plot_score_distribution(
                result.docking_poses)
        if result.metal_scores:
            plots["energy_decomp"] = self.visualizer.plot_energy_decomposition(
                result.metal_scores)
        if result.metal_scores:
            plots["radar"] = self.visualizer.plot_metal_geometry_radar(
                result.metal_scores[0])
        if result.certification:
            plots["certification"] = \
                self.visualizer.plot_certification_convergence(
                    result.certification)
        if result.rescoring:
            plots["funnel"] = self.visualizer.plot_multifidelity_funnel(
                result.rescoring)
        return plots
