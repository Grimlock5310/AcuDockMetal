"""
Visualization Module — 3D interactive views (py3Dmol) and 2D analysis
plots (matplotlib) for metal-aware docking results.
"""

from __future__ import annotations

import logging
from typing import Dict, List, Optional, Tuple

import numpy as np

log = logging.getLogger(__name__)

# Color palettes
METAL_COLORS = {
    "ZN": "#7C7C7C", "FE": "#E06633", "CU": "#C88033",
    "CO": "#F090A0", "NI": "#50D050", "MN": "#9C7AC7",
    "MG": "#8AFF00", "CA": "#3DFF00", "LN": "#70D4FF",
}

HYPOTHESIS_COLORS = [
    "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728",
    "#9467bd", "#8c564b", "#e377c2", "#7f7f7f",
]

DONOR_COLORS = {"N": "#3050F8", "O": "#FF0D0D", "S": "#FFFF30"}


class DockingVisualizer:
    """
    Visualization tools for AcuDockMetal results.

    3D views via py3Dmol (for Jupyter/Colab) and 2D plots via matplotlib.
    """

    def __init__(self, width: int = 800, height: int = 500) -> None:
        self.width = width
        self.height = height

    # ------------------------------------------------------------------
    # 3D Visualization (py3Dmol)
    # ------------------------------------------------------------------

    def view_metal_site(
        self,
        pdb_path: str,
        metal_site,
        show_waters: bool = True,
    ):
        """
        3D view of a metal site in protein context.

        Shows protein as cartoon, metal as sphere, coordination
        bonds as dashed lines, and donor residues as sticks.
        """
        import py3Dmol

        viewer = py3Dmol.view(width=self.width, height=self.height)

        with open(pdb_path) as f:
            pdb_str = f.read()
        viewer.addModel(pdb_str, "pdb")

        # Protein cartoon
        viewer.setStyle({"cartoon": {"color": "spectrum", "opacity": 0.7}})

        # Metal sphere
        metal_sel = {
            "resn": metal_site.pdb_resname,
            "resi": metal_site.residue_number,
        }
        color = METAL_COLORS.get(metal_site.metal_symbol, "#AAAAAA")
        viewer.setStyle(metal_sel, {
            "sphere": {"radius": 1.0, "color": color}
        })

        # Coordinating residues as sticks
        for atom in metal_site.coordinating_atoms:
            if not show_waters and atom.is_water:
                continue
            res_sel = {
                "chain": atom.chain_id,
                "resi": atom.residue_number,
                "resn": atom.residue_name,
            }
            viewer.setStyle(res_sel, {
                "stick": {"radius": 0.15},
                "cartoon": {"color": "spectrum", "opacity": 0.7},
            })

            # Coordination bond (dashed line)
            d_color = DONOR_COLORS.get(atom.element, "#FFFFFF")
            viewer.addLine({
                "start": {
                    "x": float(metal_site.coord[0]),
                    "y": float(metal_site.coord[1]),
                    "z": float(metal_site.coord[2]),
                },
                "end": {
                    "x": float(atom.coord[0]),
                    "y": float(atom.coord[1]),
                    "z": float(atom.coord[2]),
                },
                "color": d_color,
                "dashed": True,
            })

        # Labels
        viewer.addLabel(
            metal_site.label,
            {
                "position": {
                    "x": float(metal_site.coord[0]),
                    "y": float(metal_site.coord[1]),
                    "z": float(metal_site.coord[2] + 1.5),
                },
                "fontSize": 12,
                "fontColor": "white",
                "backgroundColor": color,
            },
        )

        viewer.zoomTo(metal_sel)
        return viewer

    def view_docked_pose(
        self,
        pdb_path: str,
        pose,
        hypothesis=None,
        ligand_pdb_block: Optional[str] = None,
    ):
        """
        3D view of a docked pose in protein context.

        Shows protein, metal, docked ligand, and coordination bonds.
        """
        import py3Dmol

        viewer = py3Dmol.view(width=self.width, height=self.height)

        with open(pdb_path) as f:
            pdb_str = f.read()
        viewer.addModel(pdb_str, "pdb")
        viewer.setStyle({"cartoon": {"color": "spectrum", "opacity": 0.5}})

        # Add ligand
        if ligand_pdb_block:
            viewer.addModel(ligand_pdb_block, "pdb")
            viewer.setStyle(
                {"model": 1},
                {"stick": {"colorscheme": "greenCarbon", "radius": 0.2}},
            )

        # Metal site if hypothesis provided
        if hypothesis is not None:
            site = hypothesis.metal_site
            color = METAL_COLORS.get(site.metal_symbol, "#AAAAAA")
            viewer.setStyle(
                {"resn": site.pdb_resname, "resi": site.residue_number},
                {"sphere": {"radius": 0.8, "color": color}},
            )

            # Show coordination lines from metal to nearby donors
            for donor in hypothesis.protein_donors:
                viewer.addLine({
                    "start": {
                        "x": float(site.coord[0]),
                        "y": float(site.coord[1]),
                        "z": float(site.coord[2]),
                    },
                    "end": {
                        "x": float(donor.coord[0]),
                        "y": float(donor.coord[1]),
                        "z": float(donor.coord[2]),
                    },
                    "color": "#AAAAAA",
                    "dashed": True,
                })

        viewer.zoomTo()
        return viewer

    def view_pose_comparison(
        self,
        pdb_path: str,
        poses: list,
        ligand_pdb_blocks: Optional[List[str]] = None,
        max_poses: int = 5,
    ):
        """Overlay multiple docked poses with distinct colors."""
        import py3Dmol

        viewer = py3Dmol.view(width=self.width, height=self.height)

        with open(pdb_path) as f:
            pdb_str = f.read()
        viewer.addModel(pdb_str, "pdb")
        viewer.setStyle({"cartoon": {"color": "spectrum", "opacity": 0.4}})

        if ligand_pdb_blocks is None:
            ligand_pdb_blocks = []

        for i, block in enumerate(ligand_pdb_blocks[:max_poses]):
            color = HYPOTHESIS_COLORS[i % len(HYPOTHESIS_COLORS)]
            viewer.addModel(block, "pdb")
            viewer.setStyle(
                {"model": i + 1},
                {"stick": {"color": color, "radius": 0.15}},
            )

        viewer.zoomTo()
        return viewer

    # ------------------------------------------------------------------
    # 2D Plots (matplotlib)
    # ------------------------------------------------------------------

    def plot_score_distribution(
        self,
        poses: list,
        title: str = "Docking Score Distribution",
    ):
        """Histogram of Vina scores colored by hypothesis."""
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(8, 4))

        hyp_ids = sorted(set(p.hypothesis_id for p in poses))
        for hi in hyp_ids:
            scores = [p.vina_score for p in poses if p.hypothesis_id == hi]
            color = HYPOTHESIS_COLORS[hi % len(HYPOTHESIS_COLORS)]
            ax.hist(scores, bins=20, alpha=0.6, color=color,
                    label=f"Hypothesis {hi}")

        ax.set_xlabel("Vina Score (kcal/mol)")
        ax.set_ylabel("Count")
        ax.set_title(title)
        ax.legend()
        plt.tight_layout()
        return fig

    def plot_metal_geometry_radar(
        self,
        score_result,
        title: str = "Metal Coordination Quality",
    ):
        """Radar plot of metal scoring components."""
        import matplotlib.pyplot as plt

        categories = ["Distance", "Angle", "CN", "Donor", "Geometry"]
        values = [
            score_result.breakdown.distance_penalty,
            score_result.breakdown.angle_penalty,
            score_result.breakdown.cn_penalty,
            score_result.breakdown.donor_compatibility,
            score_result.breakdown.geometry_rmsd,
        ]

        # Normalize to 0-1 for radar display
        max_val = max(values) if max(values) > 0 else 1.0
        values_norm = [v / max_val for v in values]
        values_norm.append(values_norm[0])  # close the polygon

        angles = np.linspace(0, 2 * np.pi, len(categories),
                             endpoint=False).tolist()
        angles.append(angles[0])

        fig, ax = plt.subplots(figsize=(6, 6),
                               subplot_kw=dict(polar=True))
        ax.fill(angles, values_norm, alpha=0.25, color="#1f77b4")
        ax.plot(angles, values_norm, "o-", color="#1f77b4")
        ax.set_xticks(angles[:-1])
        ax.set_xticklabels(categories)
        ax.set_title(title)
        plt.tight_layout()
        return fig

    def plot_certification_convergence(
        self,
        cert_result,
        title: str = "Certification Convergence",
    ):
        """Plot gap vs nodes explored during branch-and-bound."""
        import matplotlib.pyplot as plt

        if not cert_result.gap_history:
            fig, ax = plt.subplots(figsize=(8, 4))
            ax.text(0.5, 0.5, "No convergence data available",
                    ha="center", va="center", transform=ax.transAxes)
            ax.set_title(title)
            return fig

        nodes, gaps = zip(*cert_result.gap_history)

        fig, ax = plt.subplots(figsize=(8, 4))
        ax.plot(nodes, gaps, "b-", linewidth=2)
        ax.axhline(y=cert_result.epsilon, color="r", linestyle="--",
                    label=f"epsilon = {cert_result.epsilon}")
        ax.set_xlabel("Nodes Explored")
        ax.set_ylabel("Gap (kcal/mol)")
        ax.set_title(title)
        ax.legend()
        ax.set_yscale("log")
        plt.tight_layout()
        return fig

    def plot_multifidelity_funnel(
        self,
        rescore_result,
        title: str = "Multi-Fidelity Rescoring Funnel",
    ):
        """Funnel diagram showing pose filtering at each level."""
        import matplotlib.pyplot as plt

        levels = ["Input", "L1: Metal", "L2: CNN", "L3: MM", "L4: QM/MM"]
        counts = [
            rescore_result.total_input,
            rescore_result.level1_survivors,
            rescore_result.level2_survivors,
            rescore_result.level3_survivors,
            rescore_result.level4_survivors,
        ]

        fig, ax = plt.subplots(figsize=(8, 5))

        colors = ["#d4e6f1", "#85c1e9", "#3498db", "#2471a3", "#1a5276"]
        max_width = 1.0
        for i, (level, count) in enumerate(zip(levels, counts)):
            width = max_width * (count / max(counts[0], 1))
            width = max(width, 0.05)  # minimum visibility
            ax.barh(
                len(levels) - 1 - i, width,
                height=0.6, left=(max_width - width) / 2,
                color=colors[i], edgecolor="white",
            )
            ax.text(max_width / 2, len(levels) - 1 - i,
                    f"{level}\n{count} poses",
                    ha="center", va="center", fontsize=10,
                    fontweight="bold")

        ax.set_xlim(0, max_width)
        ax.set_ylim(-0.5, len(levels) - 0.5)
        ax.axis("off")
        ax.set_title(title, fontsize=14)
        plt.tight_layout()
        return fig

    def plot_energy_decomposition(
        self,
        score_results: list,
        title: str = "Score Component Decomposition",
    ):
        """Stacked bar chart of metal score components per pose."""
        import matplotlib.pyplot as plt

        if not score_results:
            fig, ax = plt.subplots()
            ax.text(0.5, 0.5, "No data", ha="center", va="center")
            return fig

        pose_labels = [f"P{r.pose_id}" for r in score_results]
        components = {
            "Distance": [r.breakdown.distance_penalty for r in score_results],
            "Angle": [r.breakdown.angle_penalty for r in score_results],
            "CN": [r.breakdown.cn_penalty for r in score_results],
            "Donor": [r.breakdown.donor_compatibility for r in score_results],
            "Geometry": [r.breakdown.geometry_rmsd for r in score_results],
        }

        fig, ax = plt.subplots(figsize=(10, 5))

        bottom = np.zeros(len(score_results))
        colors = ["#e74c3c", "#f39c12", "#3498db", "#2ecc71", "#9b59b6"]
        for (name, vals), color in zip(components.items(), colors):
            ax.bar(pose_labels, vals, bottom=bottom, label=name, color=color)
            bottom += np.array(vals)

        ax.set_xlabel("Pose")
        ax.set_ylabel("Penalty Score")
        ax.set_title(title)
        ax.legend()
        plt.tight_layout()
        return fig
