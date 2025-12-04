#!/usr/bin/env python3
"""
Principal Component Analysis for ancestry inference.

This script performs PCA on merged genetic data to visualize population
structure and project samples onto reference populations.

Features:
    - Use PLINK or EIGENSOFT smartpca for PCA
    - Project sample onto reference populations
    - Generate PCA plots (PC1 vs PC2, PC3 vs PC4, PC5 vs PC6)
    - Color-code by population/region
    - Save eigenvalues, eigenvectors, PC coordinates
    - Create static plots with matplotlib

Usage:
    python run_pca.py --input data/processed/merged \\
                      --output-dir results/pca/ \\
                      --n-pcs 10

Example:
    # Basic PCA
    python run_pca.py --input merged --output-dir results/pca/

    # With population labels
    python run_pca.py --input merged --output-dir results/pca/ \\
                      --population-file data/populations.txt
"""

import argparse
import logging
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from tqdm import tqdm

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from utils.config import (
    LOGS_DIR,
    RESULTS_DIR,
    PCA_DEFAULTS,
    POPULATION_COLORS,
    get_tool_path,
)
from utils.file_utils import (
    check_plink_files,
    ensure_dir,
    read_fam_file,
    read_population_file,
)


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


def run_plink_pca(
    input_prefix: Path,
    output_dir: Path,
    n_pcs: int = 10,
    plink_path: Optional[Path] = None,
) -> Tuple[bool, Optional[Path], Optional[Path]]:
    """
    Run PCA using PLINK.

    Args:
        input_prefix: Input PLINK file prefix
        output_dir: Output directory for results
        n_pcs: Number of principal components to compute
        plink_path: Path to PLINK executable

    Returns:
        Tuple of (success, eigenvec_path, eigenval_path)
    """
    if plink_path is None:
        plink_path = get_tool_path("plink")
        if plink_path is None:
            logger.error("PLINK not found")
            return False, None, None

    output_prefix = output_dir / "pca"

    cmd = [
        str(plink_path),
        "--bfile", str(input_prefix),
        "--pca", str(n_pcs),
        "--out", str(output_prefix),
        "--allow-no-sex",
    ]

    logger.info(f"Running PCA with PLINK ({n_pcs} components)...")

    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=False,
        )

        if result.returncode != 0:
            logger.error(f"PLINK PCA failed: {result.stderr}")
            return False, None, None

        eigenvec_path = output_prefix.parent / f"{output_prefix.name}.eigenvec"
        eigenval_path = output_prefix.parent / f"{output_prefix.name}.eigenval"

        if eigenvec_path.exists() and eigenval_path.exists():
            logger.info("PCA completed successfully")
            return True, eigenvec_path, eigenval_path
        else:
            logger.error("PCA output files not found")
            return False, None, None

    except Exception as e:
        logger.error(f"Error running PLINK PCA: {e}")
        return False, None, None


def run_smartpca(
    input_prefix: Path,
    output_dir: Path,
    n_pcs: int = 10,
    smartpca_path: Optional[Path] = None,
) -> Tuple[bool, Optional[Path], Optional[Path]]:
    """
    Run PCA using EIGENSOFT smartpca.

    Args:
        input_prefix: Input PLINK file prefix
        output_dir: Output directory for results
        n_pcs: Number of principal components to compute
        smartpca_path: Path to smartpca executable

    Returns:
        Tuple of (success, eigenvec_path, eigenval_path)
    """
    if smartpca_path is None:
        smartpca_path = get_tool_path("smartpca")
        if smartpca_path is None:
            logger.warning("smartpca not found, falling back to PLINK PCA")
            return run_plink_pca(input_prefix, output_dir, n_pcs)

    # Create parameter file
    par_file = output_dir / "smartpca.par"
    evec_file = output_dir / "pca.evec"
    eval_file = output_dir / "pca.eval"

    with open(par_file, "w") as f:
        f.write(f"genotypename: {input_prefix}.bed\n")
        f.write(f"snpname: {input_prefix}.bim\n")
        f.write(f"indivname: {input_prefix}.fam\n")
        f.write(f"evecoutname: {evec_file}\n")
        f.write(f"evaloutname: {eval_file}\n")
        f.write(f"numoutevec: {n_pcs}\n")
        f.write(f"numoutlieriter: 0\n")
        f.write(f"numthreads: 4\n")

    logger.info(f"Running PCA with smartpca ({n_pcs} components)...")

    try:
        result = subprocess.run(
            [str(smartpca_path), "-p", str(par_file)],
            capture_output=True,
            text=True,
            check=False,
        )

        if result.returncode != 0:
            logger.warning(f"smartpca failed, falling back to PLINK: {result.stderr}")
            return run_plink_pca(input_prefix, output_dir, n_pcs)

        if evec_file.exists() and eval_file.exists():
            # Convert to PLINK-like format
            eigenvec_path = output_dir / "pca.eigenvec"
            eigenval_path = output_dir / "pca.eigenval"

            convert_smartpca_output(evec_file, eval_file, eigenvec_path, eigenval_path)

            logger.info("PCA completed successfully")
            return True, eigenvec_path, eigenval_path
        else:
            logger.warning("smartpca output not found, falling back to PLINK")
            return run_plink_pca(input_prefix, output_dir, n_pcs)

    except Exception as e:
        logger.warning(f"smartpca error, falling back to PLINK: {e}")
        return run_plink_pca(input_prefix, output_dir, n_pcs)


def convert_smartpca_output(
    evec_file: Path,
    eval_file: Path,
    eigenvec_out: Path,
    eigenval_out: Path,
) -> None:
    """Convert smartpca output to PLINK-like format."""
    # Read eigenvalues
    eigenvals = []
    with open(eval_file, "r") as f:
        for line in f:
            val = line.strip()
            if val:
                eigenvals.append(float(val))

    # Write eigenvalues
    with open(eigenval_out, "w") as f:
        for val in eigenvals:
            f.write(f"{val}\n")

    # Read and convert eigenvectors
    with open(evec_file, "r") as f_in, open(eigenvec_out, "w") as f_out:
        for line in f_in:
            parts = line.strip().split()
            if len(parts) > 1 and not parts[0].startswith("#"):
                # smartpca format: IID PC1 PC2 ... population
                sample_id = parts[0]
                pcs = parts[1:-1]  # Exclude last column (population)
                f_out.write(f"{sample_id} {sample_id} {' '.join(pcs)}\n")


def load_pca_results(
    eigenvec_file: Path,
    eigenval_file: Path,
) -> Tuple[pd.DataFrame, np.ndarray]:
    """
    Load PCA results from PLINK output files.

    Args:
        eigenvec_file: Path to .eigenvec file
        eigenval_file: Path to .eigenval file

    Returns:
        Tuple of (eigenvectors DataFrame, eigenvalues array)
    """
    # Load eigenvalues
    eigenvalues = np.loadtxt(eigenval_file)

    # Load eigenvectors
    eigenvectors = pd.read_csv(
        eigenvec_file,
        sep=r"\s+",
        header=None,
    )

    # Set column names
    n_pcs = len(eigenvectors.columns) - 2
    columns = ["FID", "IID"] + [f"PC{i+1}" for i in range(n_pcs)]
    eigenvectors.columns = columns

    return eigenvectors, eigenvalues


def add_population_labels(
    eigenvectors: pd.DataFrame,
    population_file: Optional[Path] = None,
    fam_file: Optional[Path] = None,
) -> pd.DataFrame:
    """
    Add population labels to eigenvector DataFrame.

    Args:
        eigenvectors: Eigenvector DataFrame
        population_file: Optional file mapping IID to population
        fam_file: Optional FAM file (uses FID as population)

    Returns:
        DataFrame with Population column
    """
    eigenvectors = eigenvectors.copy()

    if population_file and population_file.exists():
        pop_map = read_population_file(population_file)
        eigenvectors["Population"] = eigenvectors["IID"].map(pop_map)
        eigenvectors["Population"] = eigenvectors["Population"].fillna("Unknown")
    elif fam_file and fam_file.exists():
        # Use FID as population label
        fam_df = read_fam_file(fam_file)
        pop_map = dict(zip(fam_df["IID"], fam_df["FID"]))
        eigenvectors["Population"] = eigenvectors["IID"].map(pop_map)
        eigenvectors["Population"] = eigenvectors["Population"].fillna("Unknown")
    else:
        eigenvectors["Population"] = "Sample"

    return eigenvectors


def get_population_color(population: str) -> str:
    """Get color for a population."""
    return POPULATION_COLORS.get(population, "#808080")


def create_pca_plot(
    eigenvectors: pd.DataFrame,
    eigenvalues: np.ndarray,
    pc_x: int,
    pc_y: int,
    output_file: Path,
    highlight_sample: Optional[str] = None,
    title: Optional[str] = None,
    figsize: Tuple[int, int] = (12, 10),
) -> None:
    """
    Create a PCA scatter plot.

    Args:
        eigenvectors: DataFrame with PC coordinates and Population column
        eigenvalues: Array of eigenvalues for variance explained
        pc_x: PC number for x-axis (1-based)
        pc_y: PC number for y-axis (1-based)
        output_file: Output file path
        highlight_sample: Sample ID to highlight
        title: Plot title
        figsize: Figure size
    """
    fig, ax = plt.subplots(figsize=figsize)

    # Calculate variance explained
    total_var = np.sum(eigenvalues)
    var_x = eigenvalues[pc_x - 1] / total_var * 100 if pc_x <= len(eigenvalues) else 0
    var_y = eigenvalues[pc_y - 1] / total_var * 100 if pc_y <= len(eigenvalues) else 0

    # Get unique populations
    populations = eigenvectors["Population"].unique()

    # Plot each population
    for pop in sorted(populations):
        mask = eigenvectors["Population"] == pop
        subset = eigenvectors[mask]

        color = get_population_color(pop)
        alpha = 0.7
        size = 30

        # Check if this population contains the highlighted sample
        if highlight_sample and highlight_sample in subset["IID"].values:
            # Plot other points normally
            other_mask = subset["IID"] != highlight_sample
            ax.scatter(
                subset.loc[other_mask, f"PC{pc_x}"],
                subset.loc[other_mask, f"PC{pc_y}"],
                c=color,
                label=pop,
                alpha=alpha,
                s=size,
                edgecolors="none",
            )
            # Highlight the sample
            sample_mask = subset["IID"] == highlight_sample
            ax.scatter(
                subset.loc[sample_mask, f"PC{pc_x}"],
                subset.loc[sample_mask, f"PC{pc_y}"],
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
                subset[f"PC{pc_x}"],
                subset[f"PC{pc_y}"],
                c=color,
                label=pop,
                alpha=alpha,
                s=size,
                edgecolors="none",
            )

    # Labels and title
    ax.set_xlabel(f"PC{pc_x} ({var_x:.1f}% variance explained)", fontsize=12)
    ax.set_ylabel(f"PC{pc_y} ({var_y:.1f}% variance explained)", fontsize=12)

    if title:
        ax.set_title(title, fontsize=14)
    else:
        ax.set_title(f"PCA: PC{pc_x} vs PC{pc_y}", fontsize=14)

    # Legend
    ax.legend(
        bbox_to_anchor=(1.02, 1),
        loc="upper left",
        fontsize=8,
        ncol=1 if len(populations) <= 20 else 2,
    )

    # Grid
    ax.grid(True, alpha=0.3)

    # Save
    plt.tight_layout()
    plt.savefig(output_file, dpi=150, bbox_inches="tight")
    plt.close()

    logger.info(f"Created plot: {output_file}")


def create_variance_explained_plot(
    eigenvalues: np.ndarray,
    output_file: Path,
    n_pcs: int = 10,
) -> None:
    """
    Create a variance explained plot (scree plot).

    Args:
        eigenvalues: Array of eigenvalues
        output_file: Output file path
        n_pcs: Number of PCs to show
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    # Calculate variance explained
    total_var = np.sum(eigenvalues)
    var_explained = eigenvalues / total_var * 100
    cumulative_var = np.cumsum(var_explained)

    # Limit to n_pcs
    n_show = min(n_pcs, len(eigenvalues))
    pcs = range(1, n_show + 1)

    # Individual variance
    ax1.bar(pcs, var_explained[:n_show], color="steelblue", alpha=0.7)
    ax1.set_xlabel("Principal Component", fontsize=12)
    ax1.set_ylabel("Variance Explained (%)", fontsize=12)
    ax1.set_title("Variance Explained by Each PC", fontsize=14)
    ax1.set_xticks(pcs)

    # Cumulative variance
    ax2.plot(pcs, cumulative_var[:n_show], "o-", color="steelblue", linewidth=2)
    ax2.axhline(y=80, color="red", linestyle="--", alpha=0.7, label="80% threshold")
    ax2.set_xlabel("Principal Component", fontsize=12)
    ax2.set_ylabel("Cumulative Variance Explained (%)", fontsize=12)
    ax2.set_title("Cumulative Variance Explained", fontsize=14)
    ax2.set_xticks(pcs)
    ax2.legend()
    ax2.set_ylim(0, 100)

    plt.tight_layout()
    plt.savefig(output_file, dpi=150, bbox_inches="tight")
    plt.close()

    logger.info(f"Created variance plot: {output_file}")


def save_pca_summary(
    eigenvectors: pd.DataFrame,
    eigenvalues: np.ndarray,
    output_dir: Path,
) -> None:
    """
    Save PCA summary statistics.

    Args:
        eigenvectors: DataFrame with PC coordinates
        eigenvalues: Array of eigenvalues
        output_dir: Output directory
    """
    # Save formatted eigenvalues
    total_var = np.sum(eigenvalues)
    var_explained = eigenvalues / total_var * 100
    cumulative_var = np.cumsum(var_explained)

    with open(output_dir / "eigenvalues.txt", "w") as f:
        f.write("PC\tEigenvalue\tVariance_Explained_%\tCumulative_%\n")
        for i, (ev, ve, cv) in enumerate(zip(eigenvalues, var_explained, cumulative_var)):
            f.write(f"{i+1}\t{ev:.6f}\t{ve:.4f}\t{cv:.4f}\n")

    # Save eigenvectors with nice formatting
    eigenvectors.to_csv(output_dir / "eigenvectors.txt", sep="\t", index=False)

    logger.info(f"Saved PCA summary to {output_dir}")


def run_pca_pipeline(
    input_prefix: Path,
    output_dir: Path,
    n_pcs: int = 10,
    population_file: Optional[Path] = None,
    highlight_sample: Optional[str] = None,
    use_smartpca: bool = False,
) -> bool:
    """
    Run the complete PCA pipeline.

    Args:
        input_prefix: Input PLINK file prefix
        output_dir: Output directory
        n_pcs: Number of principal components
        population_file: Optional population labels file
        highlight_sample: Optional sample ID to highlight
        use_smartpca: Whether to try smartpca first

    Returns:
        True if successful
    """
    ensure_dir(output_dir)

    # Run PCA
    if use_smartpca:
        success, eigenvec_path, eigenval_path = run_smartpca(
            input_prefix, output_dir, n_pcs
        )
    else:
        success, eigenvec_path, eigenval_path = run_plink_pca(
            input_prefix, output_dir, n_pcs
        )

    if not success or eigenvec_path is None or eigenval_path is None:
        logger.error("PCA failed")
        return False

    # Load results
    eigenvectors, eigenvalues = load_pca_results(eigenvec_path, eigenval_path)

    # Add population labels
    fam_file = input_prefix.parent / f"{input_prefix.name}.fam"
    eigenvectors = add_population_labels(eigenvectors, population_file, fam_file)

    # Create plots
    logger.info("Creating PCA plots...")

    # PC1 vs PC2
    create_pca_plot(
        eigenvectors, eigenvalues, 1, 2,
        output_dir / "pca_plot_PC1_PC2.png",
        highlight_sample=highlight_sample,
        title="PCA: PC1 vs PC2",
    )

    # PC3 vs PC4
    if n_pcs >= 4:
        create_pca_plot(
            eigenvectors, eigenvalues, 3, 4,
            output_dir / "pca_plot_PC3_PC4.png",
            highlight_sample=highlight_sample,
            title="PCA: PC3 vs PC4",
        )

    # PC5 vs PC6
    if n_pcs >= 6:
        create_pca_plot(
            eigenvectors, eigenvalues, 5, 6,
            output_dir / "pca_plot_PC5_PC6.png",
            highlight_sample=highlight_sample,
            title="PCA: PC5 vs PC6",
        )

    # Variance explained plot
    create_variance_explained_plot(
        eigenvalues,
        output_dir / "variance_explained.png",
        n_pcs=n_pcs,
    )

    # Save summary
    save_pca_summary(eigenvectors, eigenvalues, output_dir)

    # Clean up PLINK log files
    for log_file in output_dir.glob("*.log"):
        log_file.unlink(missing_ok=True)
    for nosex_file in output_dir.glob("*.nosex"):
        nosex_file.unlink(missing_ok=True)

    logger.info("PCA pipeline completed successfully")
    return True


def main() -> int:
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Principal Component Analysis for ancestry inference.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Basic PCA
    %(prog)s --input merged --output-dir results/pca/

    # With population labels
    %(prog)s --input merged --output-dir results/pca/ \\
             --population-file populations.txt

    # Highlight a specific sample
    %(prog)s --input merged --output-dir results/pca/ \\
             --sample-id MYSAMPLE

Output:
    - eigenvalues.txt - Eigenvalues with variance explained
    - eigenvectors.txt - PC coordinates for all samples
    - pca_plot_PC1_PC2.png - PC1 vs PC2 scatter plot
    - pca_plot_PC3_PC4.png - PC3 vs PC4 scatter plot
    - variance_explained.png - Scree plot
        """,
    )

    parser.add_argument(
        "--input", "-i",
        type=str,
        required=True,
        help="Input PLINK file prefix (required)",
    )

    parser.add_argument(
        "--output-dir", "-o",
        type=str,
        default=str(RESULTS_DIR / "pca"),
        help=f"Output directory (default: {RESULTS_DIR / 'pca'})",
    )

    parser.add_argument(
        "--n-pcs",
        type=int,
        default=PCA_DEFAULTS["n_pcs"],
        help=f"Number of PCs to compute (default: {PCA_DEFAULTS['n_pcs']})",
    )

    parser.add_argument(
        "--population-file",
        type=str,
        help="File mapping sample IDs to populations (optional)",
    )

    parser.add_argument(
        "--sample-id",
        type=str,
        help="Sample ID to highlight in plots (optional)",
    )

    parser.add_argument(
        "--use-smartpca",
        action="store_true",
        help="Try to use EIGENSOFT smartpca (falls back to PLINK)",
    )

    parser.add_argument(
        "--log-file",
        type=str,
        help="Path to log file",
    )

    args = parser.parse_args()

    # Set up logging
    log_file = Path(args.log_file) if args.log_file else LOGS_DIR / "run_pca.log"
    setup_logging(log_file)

    logger.info("=" * 60)
    logger.info("Principal Component Analysis")
    logger.info("=" * 60)

    # Validate input
    input_prefix = Path(args.input).resolve()
    if not check_plink_files(input_prefix):
        logger.error(f"Input PLINK files not found: {input_prefix}")
        return 1

    output_dir = Path(args.output_dir).resolve()
    population_file = Path(args.population_file) if args.population_file else None

    logger.info(f"Input: {input_prefix}")
    logger.info(f"Output directory: {output_dir}")
    logger.info(f"Number of PCs: {args.n_pcs}")

    # Run PCA pipeline
    success = run_pca_pipeline(
        input_prefix=input_prefix,
        output_dir=output_dir,
        n_pcs=args.n_pcs,
        population_file=population_file,
        highlight_sample=args.sample_id,
        use_smartpca=args.use_smartpca,
    )

    if success:
        logger.info("=" * 60)
        logger.info("PCA completed successfully!")
        logger.info(f"Results saved to: {output_dir}")
        logger.info("=" * 60)
        return 0
    else:
        logger.error("PCA failed")
        return 1


if __name__ == "__main__":
    sys.exit(main())
