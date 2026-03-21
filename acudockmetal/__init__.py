"""
AcuDockMetal: Metal-Aware epsilon-Certified Molecular Docking System

A novel docking system combining:
- Metal-aware scoring with coordination chemistry modeling
- epsilon-certified optimality guarantees via branch-and-bound
- Multi-fidelity rescoring (classical -> ML -> QM/MM)
"""

__version__ = "0.1.0"

from acudockmetal.metal_params import MetalParameterLibrary
from acudockmetal.chelator_detect import ChelatingGroupDetector
from acudockmetal.metal_site import MetalSiteDetector, CoordinationHypothesisGenerator
from acudockmetal.preparation import ReceptorPreparator, LigandPreparator
from acudockmetal.docking_engines import VinaEngine, GninaEngine, DockingOrchestrator
from acudockmetal.metal_scoring import MetalAwareScorer
from acudockmetal.certified_dock import CertifiedDockingEngine
from acudockmetal.rescoring import MultiFidelityRescorer
from acudockmetal.validation import GeometryValidator
from acudockmetal.visualization import DockingVisualizer
from acudockmetal.pipeline import AcuDockMetalPipeline
