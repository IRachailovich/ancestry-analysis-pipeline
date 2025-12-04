#!/usr/bin/env python3
"""
Create publication-quality ancestry visualization plots.

This script creates combined visualizations from PCA and ADMIXTURE results,
suitable for publication or presentation.

Features:
    - PCA plots with population labels and legends
    - ADMIXTURE bar plots with proper colors
    - Combined visualization (PCA + ADMIXTURE side-by-side)
    - Optional geographic maps if coordinates provided

Usage:
    python visualize_ancestry.py --pca-dir results/pca/ \\
                                 --admixture-dir results/admixture/ \\
                                 --output-dir results/figures/

Example:
    # Combined visualization
    python visualize_ancestry.py --pca-dir results/pca/ \\
                                 --admixture-dir results/admixture/ \\
                                 --output-dir results/figures/ \\
                                 --sample-id MYSAMPLE
"""

import argparse
import logging
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import pandas as pd

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from utils.config import (
    LOGS_DIR,
    RESULTS_DIR,
    POPULATION_COLORS,
)
from utils.file_utils import ensure_dir, read_population_file


logger = logging.getLogger(__name__)


def setup_logging(log_file: Optional[Path] = None) -> logging.Logger:
    """Set up logging configuration."""
    log_format = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"

    handlers: List[logging.Handler] = [logging.StreamHandler(sys.stdout)]

    if log_file:
        log_file.parent.mkdir(parents=True, exist_ok=True)
        handlers.append(logging.FileHandler(log_file))

    logging.basicConfig(
        level=logging.INFO,
        format=log_format,
        handlers=handlers,
    )

    return logging.getLogger(__name__)


def load_pca_results(pca_dir: Path) -> Tuple[Optional[pd.DataFrame], Optional[np.ndarray]]:
    """
    Load PCA results from directory.

    Args:
        pca_dir: PCA results directory

    Returns:
        Tuple of (eigenvectors DataFrame, eigenvalues array)
    """
    eigenvec_file = pca_dir / "eigenvectors.txt"
    eigenval_file = pca_dir / "eigenvalues.txt"

    if not eigenvec_file.exists():
        # Try PLINK format
        eigenvec_file = pca_dir / "pca.eigenvec"
        eigenval_file = pca_dir / "pca.eigenval"

    if not eigenvec_file.exists():
        logger.error(f"PCA eigenvector file not found in {pca_dir}")
        return None, None

    # Load eigenvectors
    eigenvectors = pd.read_csv(eigenvec_file, sep="\t")

    # Load eigenvalues
    if eigenval_file.exists():
        # Check if it has header
        with open(eigenval_file, "r") as f:
            first_line = f.readline().strip()
            if "Eigenvalue" in first_line or "PC" in first_line:
                eigenval_df = pd.read_csv(eigenval_file, sep="\t")
                eigenvalues = eigenval_df.iloc[:, 1].values if eigenval_df.shape[1] > 1 else eigenval_df.iloc[:, 0].values
            else:
                eigenvalues = np.loadtxt(eigenval_file)
    else:
        eigenvalues = None

    return eigenvectors, eigenvalues


def load_admixture_results(
    admixture_dir: Path,
    k: Optional[int] = None,
) -> Tuple[Optional[pd.DataFrame], Optional[int]]:
    """
    Load ADMIXTURE results from directory.

    Args:
        admixture_dir: ADMIXTURE results directory
        k: Specific K value to load (if None, uses optimal K from CV)

    Returns:
        Tuple of (Q matrix DataFrame, K value)
    """
    if k is None:
        # Find optimal K from CV errors
        cv_file = admixture_dir / "cv_errors.txt"
        if cv_file.exists():
            cv_df = pd.read_csv(cv_file, sep="\t")
            k = int(cv_df.loc[cv_df["CV_Error"].idxmin(), "K"])
            logger.info(f"Using optimal K={k} from cross-validation")
        else:
            # Find any Q file
            q_files = list(admixture_dir.glob("K*.Q"))
            if q_files:
                k = int(q_files[0].stem[1:])
            else:
                logger.error("No ADMIXTURE Q files found")
                return None, None

    q_file = admixture_dir / f"K{k}.Q"
    if not q_file.exists():
        logger.error(f"ADMIXTURE Q file not found: {q_file}")
        return None, None

    q_df = pd.read_csv(q_file, sep=r"\s+", header=None)
    return q_df, k


def get_population_color(population: str) -> str:
    """Get color for a population."""
    return POPULATION_COLORS.get(population, "#808080")


def create_pca_subplot(
    ax: plt.Axes,
    eigenvectors: pd.DataFrame,
    eigenvalues: Optional[np.ndarray],
    pc_x: int,
    pc_y: int,
    highlight_sample: Optional[str] = None,
) -> None:
    """
    Create a PCA scatter plot on a subplot.

    Args:
        ax: Matplotlib axes
        eigenvectors: DataFrame with PC coordinates
        eigenvalues: Array of eigenvalues
        pc_x: PC number for x-axis (1-based)
        pc_y: PC number for y-axis (1-based)
        highlight_sample: Sample ID to highlight
    """
    # Get population column
    pop_col = "Population" if "Population" in eigenvectors.columns else None

    if pop_col:
        populations = eigenvectors[pop_col].unique()
    else:
        populations = ["Sample"]
        eigenvectors = eigenvectors.copy()
        eigenvectors["Population"] = "Sample"
        pop_col = "Population"

    # Plot each population
    for pop in sorted(populations):
        mask = eigenvectors[pop_col] == pop
        subset = eigenvectors[mask]

        color = get_population_color(pop)
        alpha = 0.7
        size = 30

        pc_x_col = f"PC{pc_x}"
        pc_y_col = f"PC{pc_y}"

        if pc_x_col not in subset.columns or pc_y_col not in subset.columns:
            continue

        # Check if this population contains the highlighted sample
        if highlight_sample and "IID" in subset.columns and highlight_sample in subset["IID"].values:
            other_mask = subset["IID"] != highlight_sample
            ax.scatter(
                subset.loc[other_mask, pc_x_col],
                subset.loc[other_mask, pc_y_col],
                c=color,
                label=pop,
                alpha=alpha,
                s=size,
                edgecolors="none",
            )
            sample_mask = subset["IID"] == highlight_sample
            ax.scatter(
                subset.loc[sample_mask, pc_x_col],
                subset.loc[sample_mask, pc_y_col],
                c=color,
                s=200,
                marker="*",
                edgecolors="black",
                linewidths=2,
                zorder=10,
                label=f"{highlight_sample} (You)",
            )
        else:
            ax.scatter(
                subset[pc_x_col],
                subset[pc_y_col],
                c=color,
                label=pop,
                alpha=alpha,
                s=size,
                edgecolors="none",
            )

    # Calculate variance explained
    if eigenvalues is not None and len(eigenvalues) >= max(pc_x, pc_y):
        total_var = np.sum(eigenvalues)
        var_x = eigenvalues[pc_x - 1] / total_var * 100
        var_y = eigenvalues[pc_y - 1] / total_var * 100
        ax.set_xlabel(f"PC{pc_x} ({var_x:.1f}%)", fontsize=10)
        ax.set_ylabel(f"PC{pc_y} ({var_y:.1f}%)", fontsize=10)
    else:
        ax.set_xlabel(f"PC{pc_x}", fontsize=10)
        ax.set_ylabel(f"PC{pc_y}", fontsize=10)

    ax.set_title(f"PC{pc_x} vs PC{pc_y}", fontsize=12)
    ax.grid(True, alpha=0.3)


def create_admixture_subplot(
    ax: plt.Axes,
    q_df: pd.DataFrame,
    k: int,
    sample_ids: Optional[List[str]] = None,
    populations: Optional[List[str]] = None,
    highlight_sample: Optional[str] = None,
) -> None:
    """
    Create an ADMIXTURE bar plot on a subplot.

    Args:
        ax: Matplotlib axes
        q_df: DataFrame with ancestry proportions
        k: K value
        sample_ids: List of sample IDs
        populations: List of population labels
        highlight_sample: Sample ID to highlight
    """
    n_samples = len(q_df)

    # Sort by population if available
    if populations:
        q_df = q_df.copy()
        q_df["Population"] = populations
        q_df = q_df.sort_values("Population")
        sample_order = q_df.index.tolist()
        sorted_pops = q_df["Population"].values
    else:
        sample_order = list(range(n_samples))
        sorted_pops = None

    # Define colors
    component_colors = [
        "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
        "#FFFF33", "#A65628", "#F781BF", "#999999", "#66C2A5",
        "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F",
    ]

    # Plot stacked bars
    x = np.arange(n_samples)
    bottom = np.zeros(n_samples)

    for i in range(k):
        if i >= len(q_df.columns):
            break
        values = q_df.iloc[sample_order, i].values
        ax.bar(
            x, values, bottom=bottom,
            color=component_colors[i % len(component_colors)],
            width=1.0,
        )
        bottom += values

    # Highlight sample if specified
    if highlight_sample and sample_ids:
        try:
            sample_idx = sample_ids.index(highlight_sample)
            # Find position in sorted order
            if sample_order:
                pos = sample_order.index(sample_idx)
                ax.axvline(x=pos, color="black", linewidth=2, linestyle="--")
        except (ValueError, IndexError):
            pass

    ax.set_xlim(-0.5, n_samples - 0.5)
    ax.set_ylim(0, 1)
    ax.set_ylabel("Proportion", fontsize=10)
    ax.set_title(f"ADMIXTURE K={k}", fontsize=12)
    ax.set_xticks([])


def create_combined_figure(
    eigenvectors: pd.DataFrame,
    eigenvalues: Optional[np.ndarray],
    q_df: pd.DataFrame,
    k: int,
    output_file: Path,
    sample_ids: Optional[List[str]] = None,
    populations: Optional[List[str]] = None,
    highlight_sample: Optional[str] = None,
) -> None:
    """
    Create a combined PCA + ADMIXTURE figure.

    Args:
        eigenvectors: DataFrame with PC coordinates
        eigenvalues: Array of eigenvalues
        q_df: DataFrame with ancestry proportions
        k: ADMIXTURE K value
        output_file: Output file path
        sample_ids: List of sample IDs
        populations: List of population labels
        highlight_sample: Sample ID to highlight
    """
    fig = plt.figure(figsize=(18, 12))

    # Create grid: 2x3 layout
    # Top row: PCA plots (PC1vPC2, PC3vPC4, legend)
    # Bottom row: ADMIXTURE bar plot spanning all columns

    gs = gridspec.GridSpec(2, 3, height_ratios=[1, 0.5], hspace=0.3, wspace=0.3)

    # PCA plot 1: PC1 vs PC2
    ax1 = fig.add_subplot(gs[0, 0])
    create_pca_subplot(ax1, eigenvectors, eigenvalues, 1, 2, highlight_sample)

    # PCA plot 2: PC3 vs PC4
    ax2 = fig.add_subplot(gs[0, 1])
    create_pca_subplot(ax2, eigenvectors, eigenvalues, 3, 4, highlight_sample)

    # Legend
    ax3 = fig.add_subplot(gs[0, 2])
    ax3.axis("off")

    # Get unique populations for legend
    if "Population" in eigenvectors.columns:
        pops = sorted(eigenvectors["Population"].unique())
        handles = []
        labels = []
        for pop in pops:
            color = get_population_color(pop)
            handles.append(plt.scatter([], [], c=color, s=50))
            labels.append(pop)

        ax3.legend(
            handles, labels,
            loc="center",
            fontsize=8,
            ncol=2 if len(pops) > 15 else 1,
            title="Populations",
        )

    # ADMIXTURE bar plot spanning bottom row
    ax4 = fig.add_subplot(gs[1, :])
    create_admixture_subplot(ax4, q_df, k, sample_ids, populations, highlight_sample)

    # Overall title
    fig.suptitle("Ancestry Analysis Summary", fontsize=16, fontweight="bold", y=0.98)

    plt.savefig(output_file, dpi=150, bbox_inches="tight")
    plt.close()

    logger.info(f"Combined figure saved to: {output_file}")


def create_individual_ancestry_report(
    eigenvectors: pd.DataFrame,
    q_df: pd.DataFrame,
    k: int,
    sample_id: str,
    output_file: Path,
) -> None:
    """
    Create an individual ancestry report for a specific sample.

    Args:
        eigenvectors: DataFrame with PC coordinates
        q_df: DataFrame with ancestry proportions
        k: ADMIXTURE K value
        sample_id: Sample ID to report on
        output_file: Output file path
    """
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # Find sample index
    if "IID" in eigenvectors.columns:
        sample_mask = eigenvectors["IID"] == sample_id
        if not sample_mask.any():
            logger.warning(f"Sample {sample_id} not found in PCA results")
            sample_idx = 0
        else:
            sample_idx = eigenvectors[sample_mask].index[0]
    else:
        sample_idx = 0

    # Pie chart of ancestry proportions
    ax1 = axes[0]
    ancestry_values = q_df.iloc[sample_idx].values[:k]

    component_colors = [
        "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
        "#FFFF33", "#A65628", "#F781BF", "#999999", "#66C2A5",
    ]

    labels = [f"Ancestry {i+1}\n({v*100:.1f}%)" for i, v in enumerate(ancestry_values)]
    colors = [component_colors[i % len(component_colors)] for i in range(k)]

    # Only show labels for significant components
    labels_display = [l if v > 0.05 else "" for l, v in zip(labels, ancestry_values)]

    ax1.pie(
        ancestry_values,
        labels=labels_display,
        colors=colors,
        autopct=lambda p: f"{p:.1f}%" if p > 5 else "",
        startangle=90,
    )
    ax1.set_title(f"Ancestry Proportions (K={k})", fontsize=12)

    # Bar chart of ancestry proportions
    ax2 = axes[1]
    x = range(k)
    ax2.bar(x, ancestry_values * 100, color=colors)
    ax2.set_xlabel("Ancestry Component", fontsize=10)
    ax2.set_ylabel("Proportion (%)", fontsize=10)
    ax2.set_title("Ancestry Breakdown", fontsize=12)
    ax2.set_xticks(x)
    ax2.set_xticklabels([f"Anc {i+1}" for i in range(k)])
    ax2.set_ylim(0, 100)

    # Add value labels on bars
    for i, v in enumerate(ancestry_values):
        ax2.text(i, v * 100 + 2, f"{v*100:.1f}%", ha="center", fontsize=9)

    fig.suptitle(f"Ancestry Report: {sample_id}", fontsize=14, fontweight="bold")

    plt.tight_layout()
    plt.savefig(output_file, dpi=150, bbox_inches="tight")
    plt.close()

    logger.info(f"Individual report saved to: {output_file}")


def run_visualization_pipeline(
    pca_dir: Optional[Path],
    admixture_dir: Optional[Path],
    output_dir: Path,
    sample_id: Optional[str] = None,
    k: Optional[int] = None,
) -> bool:
    """
    Run the complete visualization pipeline.

    Args:
        pca_dir: PCA results directory
        admixture_dir: ADMIXTURE results directory
        output_dir: Output directory for figures
        sample_id: Optional sample ID to highlight
        k: Optional K value for ADMIXTURE

    Returns:
        True if successful
    """
    ensure_dir(output_dir)

    # Load PCA results
    eigenvectors = None
    eigenvalues = None
    if pca_dir and pca_dir.exists():
        eigenvectors, eigenvalues = load_pca_results(pca_dir)
        if eigenvectors is not None:
            logger.info(f"Loaded PCA results: {len(eigenvectors)} samples")

    # Load ADMIXTURE results
    q_df = None
    if admixture_dir and admixture_dir.exists():
        q_df, k = load_admixture_results(admixture_dir, k)
        if q_df is not None:
            logger.info(f"Loaded ADMIXTURE results: K={k}, {len(q_df)} samples")

    if eigenvectors is None and q_df is None:
        logger.error("No valid results to visualize")
        return False

    # Get sample IDs and populations
    sample_ids = None
    populations = None

    if eigenvectors is not None and "IID" in eigenvectors.columns:
        sample_ids = eigenvectors["IID"].tolist()
    if eigenvectors is not None and "Population" in eigenvectors.columns:
        populations = eigenvectors["Population"].tolist()

    # Create combined figure if both results available
    if eigenvectors is not None and q_df is not None:
        create_combined_figure(
            eigenvectors, eigenvalues, q_df, k,
            output_dir / "combined_ancestry.png",
            sample_ids, populations, sample_id,
        )

    # Create individual report if sample ID specified
    if sample_id and eigenvectors is not None and q_df is not None:
        create_individual_ancestry_report(
            eigenvectors, q_df, k, sample_id,
            output_dir / f"individual_report_{sample_id}.png",
        )

    logger.info("Visualization pipeline completed")
    return True


def main() -> int:
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Create publication-quality ancestry visualization plots.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Combined visualization
    %(prog)s --pca-dir results/pca/ --admixture-dir results/admixture/ \\
             --output-dir results/figures/

    # With sample highlight
    %(prog)s --pca-dir results/pca/ --admixture-dir results/admixture/ \\
             --output-dir results/figures/ --sample-id MYSAMPLE

Output:
    - combined_ancestry.png - Combined PCA + ADMIXTURE visualization
    - individual_report_{sample_id}.png - Individual ancestry report
        """,
    )

    parser.add_argument(
        "--pca-dir",
        type=str,
        help="PCA results directory",
    )

    parser.add_argument(
        "--admixture-dir",
        type=str,
        help="ADMIXTURE results directory",
    )

    parser.add_argument(
        "--output-dir", "-o",
        type=str,
        default=str(RESULTS_DIR / "figures"),
        help=f"Output directory (default: {RESULTS_DIR / 'figures'})",
    )

    parser.add_argument(
        "--sample-id",
        type=str,
        help="Sample ID to highlight (optional)",
    )

    parser.add_argument(
        "--k",
        type=int,
        help="ADMIXTURE K value to use (default: optimal from CV)",
    )

    parser.add_argument(
        "--log-file",
        type=str,
        help="Path to log file",
    )

    args = parser.parse_args()

    # Set up logging
    log_file = Path(args.log_file) if args.log_file else LOGS_DIR / "visualize_ancestry.log"
    setup_logging(log_file)

    logger.info("=" * 60)
    logger.info("Ancestry Visualization")
    logger.info("=" * 60)

    # Validate inputs
    pca_dir = Path(args.pca_dir) if args.pca_dir else None
    admixture_dir = Path(args.admixture_dir) if args.admixture_dir else None
    output_dir = Path(args.output_dir).resolve()

    if pca_dir is None and admixture_dir is None:
        logger.error("At least one of --pca-dir or --admixture-dir is required")
        return 1

    # Run visualization
    success = run_visualization_pipeline(
        pca_dir=pca_dir,
        admixture_dir=admixture_dir,
        output_dir=output_dir,
        sample_id=args.sample_id,
        k=args.k,
    )

    if success:
        logger.info("=" * 60)
        logger.info("Visualization completed successfully!")
        logger.info(f"Figures saved to: {output_dir}")
        logger.info("=" * 60)
        return 0
    else:
        logger.error("Visualization failed")
        return 1


if __name__ == "__main__":
    sys.exit(main())
