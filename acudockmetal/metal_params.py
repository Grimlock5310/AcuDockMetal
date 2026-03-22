"""
Metal Parameter Library — comprehensive coordination chemistry data for
biologically relevant metals (transition, main-group, and lanthanide).

Data sourced from structural surveys of metalloprotein coordination
(MetalPDB, CSD), 12-6-4 LJ parameterizations, and coordination chemistry
literature as described in the AcuDockMetal research reports.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple


# ---------------------------------------------------------------------------
# Geometry templates — ideal angles for common coordination polyhedra
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class GeometryTemplate:
    """Ideal coordination geometry with reference angles."""
    name: str
    coordination_number: int
    ideal_angles: Tuple[float, ...]  # unique cis/trans angles in degrees
    description: str = ""

    def angle_deviation(self, observed_angles: List[float]) -> float:
        """Mean minimum angular deviation from ideal angles."""
        if not observed_angles:
            return 0.0
        total = 0.0
        for obs in observed_angles:
            min_dev = min(abs(obs - ideal) for ideal in self.ideal_angles)
            total += min_dev
        return total / len(observed_angles)


GEOMETRY_TEMPLATES: Dict[str, GeometryTemplate] = {
    "linear": GeometryTemplate(
        "linear", 2, (180.0,), "Linear coordination (CN=2)"),
    "trigonal_planar": GeometryTemplate(
        "trigonal_planar", 3, (120.0,), "Trigonal planar (CN=3)"),
    "tetrahedral": GeometryTemplate(
        "tetrahedral", 4, (109.5,), "Tetrahedral (CN=4)"),
    "square_planar": GeometryTemplate(
        "square_planar", 4, (90.0, 180.0), "Square planar (CN=4)"),
    "trigonal_bipyramidal": GeometryTemplate(
        "trigonal_bipyramidal", 5, (90.0, 120.0, 180.0),
        "Trigonal bipyramidal (CN=5)"),
    "square_pyramidal": GeometryTemplate(
        "square_pyramidal", 5, (90.0, 104.0, 180.0),
        "Square pyramidal (CN=5)"),
    "octahedral": GeometryTemplate(
        "octahedral", 6, (90.0, 180.0), "Octahedral (CN=6)"),
    "trigonal_prismatic": GeometryTemplate(
        "trigonal_prismatic", 6, (82.0, 136.0),
        "Trigonal prismatic (CN=6)"),
    "pentagonal_bipyramidal": GeometryTemplate(
        "pentagonal_bipyramidal", 7, (72.0, 90.0, 144.0, 180.0),
        "Pentagonal bipyramidal (CN=7)"),
    "square_antiprismatic": GeometryTemplate(
        "square_antiprismatic", 8, (70.5, 109.5, 141.6, 180.0),
        "Square antiprismatic (CN=8)"),
    "tricapped_trigonal_prismatic": GeometryTemplate(
        "tricapped_trigonal_prismatic", 9,
        (70.0, 82.0, 120.0, 136.0, 180.0),
        "Tricapped trigonal prismatic (CN=9)"),
}


# ---------------------------------------------------------------------------
# Per-metal parameter dataclass
# ---------------------------------------------------------------------------

@dataclass
class MetalParams:
    """Coordination chemistry parameters for a single metal ion type."""
    symbol: str
    name: str
    common_oxidation_states: List[int]
    preferred_geometries: Dict[int, List[str]]  # CN -> geometry names
    ideal_distances: Dict[str, float]  # donor element -> distance (angstrom)
    preferred_donors: List[str]  # ordered by preference
    typical_cn_range: Tuple[int, int]
    water_propensity: float  # 0-1: probability of coordinated waters
    lj_r_min: float  # LJ r_min/2 in angstrom (12-6-4 params)
    lj_epsilon: float  # LJ well depth in kcal/mol
    c4_coefficient: float  # 12-6-4 C4 ion-dipole term (kcal/mol * A^4)
    notes: str = ""

    @property
    def default_cn(self) -> int:
        """Most common coordination number (lowest in typical range)."""
        return self.typical_cn_range[0]

    def get_ideal_distance(self, donor_element: str) -> float:
        """Get ideal metal-donor distance; fall back to generic if unknown."""
        donor = donor_element.upper()
        if donor in self.ideal_distances:
            return self.ideal_distances[donor]
        # Generic fallback based on ionic radius sum
        return self.ideal_distances.get("O", 2.1)

    def get_geometries_for_cn(self, cn: int) -> List[GeometryTemplate]:
        """Return geometry templates for a given coordination number."""
        names = self.preferred_geometries.get(cn, [])
        return [GEOMETRY_TEMPLATES[n] for n in names if n in GEOMETRY_TEMPLATES]


# ---------------------------------------------------------------------------
# The library itself
# ---------------------------------------------------------------------------

class MetalParameterLibrary:
    """
    Central registry of metal-ion coordination parameters.

    Covers transition metals (Zn, Fe, Cu, Co, Ni, Mn), main-group cations
    (Mg, Ca), and lanthanides (generic Ln^3+).  Data is drawn from
    structural surveys (MetalPDB / CSD), AMBER 12-6-4 parameter sets,
    and metalloprotein coordination literature.
    """

    def __init__(self) -> None:
        self._metals: Dict[str, MetalParams] = {}
        self._build_library()

    # ------------------------------------------------------------------
    def _build_library(self) -> None:
        self._add_zinc()
        self._add_iron()
        self._add_copper()
        self._add_cobalt()
        self._add_nickel()
        self._add_manganese()
        self._add_magnesium()
        self._add_calcium()
        self._add_sodium()
        self._add_potassium()
        self._add_lanthanide()

    # ------------------------------------------------------------------
    # Transition metals
    # ------------------------------------------------------------------
    def _add_zinc(self) -> None:
        self._metals["ZN"] = MetalParams(
            symbol="ZN", name="Zinc",
            common_oxidation_states=[2],
            preferred_geometries={
                4: ["tetrahedral"],
                5: ["trigonal_bipyramidal", "square_pyramidal"],
                6: ["octahedral"],
            },
            ideal_distances={"N": 2.05, "O": 2.10, "S": 2.30, "CL": 2.25},
            preferred_donors=["N", "O", "S"],
            typical_cn_range=(4, 6),
            water_propensity=0.5,
            lj_r_min=1.395, lj_epsilon=0.0125,
            c4_coefficient=88.0,
            notes="Most common catalytic metal in drug targets (MMPs, HDACs, ACE).",
        )

    def _add_iron(self) -> None:
        self._metals["FE"] = MetalParams(
            symbol="FE", name="Iron",
            common_oxidation_states=[2, 3],
            preferred_geometries={
                4: ["tetrahedral"],
                5: ["square_pyramidal", "trigonal_bipyramidal"],
                6: ["octahedral"],
            },
            ideal_distances={"N": 2.10, "O": 2.05, "S": 2.25, "NE2": 2.15},
            preferred_donors=["N", "S", "O"],
            typical_cn_range=(4, 6),
            water_propensity=0.4,
            lj_r_min=1.409, lj_epsilon=0.0130,
            c4_coefficient=100.0,
            notes="Redox-active; heme and non-heme sites. QM/MM critical.",
        )

    def _add_copper(self) -> None:
        self._metals["CU"] = MetalParams(
            symbol="CU", name="Copper",
            common_oxidation_states=[1, 2],
            preferred_geometries={
                4: ["square_planar", "tetrahedral"],
                5: ["square_pyramidal", "trigonal_bipyramidal"],
                6: ["octahedral"],
            },
            ideal_distances={"N": 2.00, "O": 1.97, "S": 2.15, "NE2": 2.02},
            preferred_donors=["N", "S", "O"],
            typical_cn_range=(4, 6),
            water_propensity=0.3,
            lj_r_min=1.380, lj_epsilon=0.0120,
            c4_coefficient=92.0,
            notes="Square-planar preferred for Cu(II); Jahn-Teller distortions common.",
        )

    def _add_cobalt(self) -> None:
        self._metals["CO"] = MetalParams(
            symbol="CO", name="Cobalt",
            common_oxidation_states=[2, 3],
            preferred_geometries={
                4: ["tetrahedral"],
                5: ["square_pyramidal"],
                6: ["octahedral"],
            },
            ideal_distances={"N": 2.10, "O": 2.08, "S": 2.25},
            preferred_donors=["N", "O"],
            typical_cn_range=(4, 6),
            water_propensity=0.5,
            lj_r_min=1.404, lj_epsilon=0.0128,
            c4_coefficient=95.0,
            notes="Octahedral preference for Co(II); sometimes substitutes for Zn.",
        )

    def _add_nickel(self) -> None:
        self._metals["NI"] = MetalParams(
            symbol="NI", name="Nickel",
            common_oxidation_states=[2],
            preferred_geometries={
                4: ["square_planar", "tetrahedral"],
                5: ["square_pyramidal"],
                6: ["octahedral"],
            },
            ideal_distances={"N": 2.08, "O": 2.06, "S": 2.18},
            preferred_donors=["N", "O", "S"],
            typical_cn_range=(4, 6),
            water_propensity=0.5,
            lj_r_min=1.373, lj_epsilon=0.0115,
            c4_coefficient=86.0,
            notes="Octahedral preference for Ni(II).",
        )

    def _add_manganese(self) -> None:
        self._metals["MN"] = MetalParams(
            symbol="MN", name="Manganese",
            common_oxidation_states=[2, 3, 4],
            preferred_geometries={
                5: ["trigonal_bipyramidal", "square_pyramidal"],
                6: ["octahedral"],
            },
            ideal_distances={"N": 2.20, "O": 2.15, "S": 2.40},
            preferred_donors=["O", "N"],
            typical_cn_range=(5, 6),
            water_propensity=0.6,
            lj_r_min=1.468, lj_epsilon=0.0150,
            c4_coefficient=105.0,
            notes="O-donor preference; common in OEC, SOD.",
        )

    # ------------------------------------------------------------------
    # Main-group hard cations
    # ------------------------------------------------------------------
    def _add_magnesium(self) -> None:
        self._metals["MG"] = MetalParams(
            symbol="MG", name="Magnesium",
            common_oxidation_states=[2],
            preferred_geometries={
                5: ["trigonal_bipyramidal", "square_pyramidal"],
                6: ["octahedral"],
            },
            ideal_distances={"O": 2.07, "N": 2.15, "S": 2.55},
            preferred_donors=["O"],
            typical_cn_range=(5, 6),
            water_propensity=0.7,
            lj_r_min=1.395, lj_epsilon=0.0107,
            c4_coefficient=72.0,
            notes="Hard Lewis acid; strongly prefers O donors; frequently fully hydrated.",
        )

    def _add_calcium(self) -> None:
        self._metals["CA"] = MetalParams(
            symbol="CA", name="Calcium",
            common_oxidation_states=[2],
            preferred_geometries={
                6: ["octahedral"],
                7: ["pentagonal_bipyramidal"],
                8: ["square_antiprismatic"],
            },
            ideal_distances={"O": 2.35, "N": 2.50, "S": 2.80},
            preferred_donors=["O"],
            typical_cn_range=(6, 8),
            water_propensity=0.6,
            lj_r_min=1.713, lj_epsilon=0.0300,
            c4_coefficient=120.0,
            notes="Larger ionic radius; higher CN (6-8); O-donor dominated.",
        )

    def _add_sodium(self) -> None:
        self._metals["NA"] = MetalParams(
            symbol="NA", name="Sodium",
            common_oxidation_states=[1],
            preferred_geometries={
                5: ["trigonal_bipyramidal", "square_pyramidal"],
                6: ["octahedral"],
            },
            ideal_distances={"O": 2.40, "N": 2.50},
            preferred_donors=["O"],
            typical_cn_range=(5, 6),
            water_propensity=0.9,
            lj_r_min=1.369, lj_epsilon=0.0874,
            c4_coefficient=30.0,
            notes="Hard Lewis acid; prefers O donors; often fully hydrated.",
        )

    def _add_potassium(self) -> None:
        self._metals["K"] = MetalParams(
            symbol="K", name="Potassium",
            common_oxidation_states=[1],
            preferred_geometries={
                6: ["octahedral"],
                7: ["pentagonal_bipyramidal"],
                8: ["square_antiprismatic"],
            },
            ideal_distances={"O": 2.80, "N": 2.90},
            preferred_donors=["O"],
            typical_cn_range=(6, 8),
            water_propensity=0.9,
            lj_r_min=1.705, lj_epsilon=0.1936,
            c4_coefficient=40.0,
            notes="Large ionic radius; high CN; O-donor dominated.",
        )

    # ------------------------------------------------------------------
    # Lanthanides (generic Ln^3+)
    # ------------------------------------------------------------------
    def _add_lanthanide(self) -> None:
        self._metals["LN"] = MetalParams(
            symbol="LN", name="Lanthanide (generic)",
            common_oxidation_states=[3],
            preferred_geometries={
                8: ["square_antiprismatic"],
                9: ["tricapped_trigonal_prismatic"],
            },
            ideal_distances={"O": 2.40, "N": 2.55, "S": 2.90},
            preferred_donors=["O", "N"],
            typical_cn_range=(8, 10),
            water_propensity=0.8,
            lj_r_min=1.750, lj_epsilon=0.0350,
            c4_coefficient=150.0,
            notes="High CN (8-10); hard Lewis acid; high water retention.",
        )

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------
    def get(self, symbol: str) -> Optional[MetalParams]:
        """Look up parameters by element symbol (case-insensitive)."""
        return self._metals.get(symbol.upper())

    def __getitem__(self, symbol: str) -> MetalParams:
        key = symbol.upper()
        if key not in self._metals:
            raise KeyError(f"Metal '{symbol}' not in library.  "
                           f"Available: {list(self._metals.keys())}")
        return self._metals[key]

    def __contains__(self, symbol: str) -> bool:
        return symbol.upper() in self._metals

    @property
    def supported_metals(self) -> List[str]:
        return list(self._metals.keys())

    def get_template(self, name: str) -> Optional[GeometryTemplate]:
        """Look up a geometry template by name."""
        return GEOMETRY_TEMPLATES.get(name)

    @property
    def all_templates(self) -> Dict[str, GeometryTemplate]:
        return dict(GEOMETRY_TEMPLATES)

    def summary(self) -> str:
        """Human-readable summary table."""
        lines = [f"{'Metal':<6} {'Name':<12} {'OxSt':>5} {'CN range':>10} "
                 f"{'Preferred geom':>25} {'Donors':>12}"]
        lines.append("-" * 78)
        for m in self._metals.values():
            ox = ",".join(str(x) for x in m.common_oxidation_states)
            cn = f"{m.typical_cn_range[0]}-{m.typical_cn_range[1]}"
            geom = "; ".join(
                g for gs in m.preferred_geometries.values() for g in gs
            )[:25]
            donors = ",".join(m.preferred_donors)
            lines.append(
                f"{m.symbol:<6} {m.name:<12} {ox:>5} {cn:>10} "
                f"{geom:>25} {donors:>12}"
            )
        return "\n".join(lines)
