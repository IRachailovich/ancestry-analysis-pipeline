#!/usr/bin/env python3
"""
ADMIXTURE analysis for ancestry proportions.

This script runs ADMIXTURE for multiple K values to estimate ancestry
proportions and identifies the optimal K using cross-validation.

Features:
    - Run ADMIXTURE for K=2 to K=15 (or user-specified range)
    - Cross-validation to determine optimal K
    - Generate ancestry proportion bar plots
    - Create summary tables with ancestry proportions
    - Parallel processing for multiple K values

Usage:
    python run_admixture.py --input data/processed/merged \\
                            --output-dir results/admixture/ \\
                            --min-k 2 --max-k 10

Example:
    # Basic ADMIXTURE analysis
    python run_admixture.py --input merged --output-dir results/admixture/

    # With custom K range
    python run_admixture.py --input merged --output-dir results/admixture/ \\
                            --min-k 2 --max-k 8 --threads 8
"""

import argparse
import logging
import os
import re
import shutil
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
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
    ADMIXTURE_DEFAULTS,
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


def get_admixture_executable() -> Optional[Path]:
    """Get ADMIXTURE executable path."""
    admixture_path = get_tool_path("admixture")
    if admixture_path is None:
        logger.error("ADMIXTURE not found. Please install ADMIXTURE and ensure it's in PATH.")
        return None
    return admixture_path


def run_admixture_k(
    bed_file: Path,
    k: int,
    cv: int = 5,
    threads: int = 1,
    admixture_path: Optional[Path] = None,
    working_dir: Optional[Path] = None,
) -> Tuple[int, float, bool]:
    """
    Run ADMIXTURE for a single K value.

    Args:
        bed_file: Path to .bed file
        k: Number of ancestral populations
        cv: Cross-validation folds
        threads: Number of threads
        admixture_path: Path to ADMIXTURE executable
        working_dir: Working directory for output

    Returns:
        Tuple of (K, CV error, success)
    """
    if admixture_path is None:
        admixture_path = get_admixture_executable()
        if admixture_path is None:
            return k, float("inf"), False

    # Run from working directory to control output location
    if working_dir is None:
        working_dir = bed_file.parent

    # Build command
    cmd = [
        str(admixture_path),
        str(bed_file),
        str(k),
        f"--cv={cv}",
        f"-j{threads}",
    ]

    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=False,
            cwd=str(working_dir),
        )

        # Parse CV error from output
        cv_error = float("inf")
        cv_pattern = r"CV error \(K=\d+\): ([\d.]+)"

        for line in result.stdout.split("\n"):
            match = re.search(cv_pattern, line)
            if match:
                cv_error = float(match.group(1))
                break

        # Also check stderr
        for line in result.stderr.split("\n"):
            match = re.search(cv_pattern, line)
            if match:
                cv_error = float(match.group(1))
                break

        success = result.returncode == 0

        return k, cv_error, success

    except Exception as e:
        logger.error(f"Error running ADMIXTURE K={k}: {e}")
        return k, float("inf"), False


def run_admixture_range(
    input_prefix: Path,
    output_dir: Path,
    min_k: int = 2,
    max_k: int = 15,
    cv: int = 5,
    threads: int = 4,
) -> Tuple[Dict[int, float], int]:
    """
    Run ADMIXTURE for a range of K values.

    Args:
        input_prefix: Input PLINK file prefix
        output_dir: Output directory
        min_k: Minimum K value
        max_k: Maximum K value
        cv: Cross-validation folds
        threads: Number of threads per run

    Returns:
        Tuple of (cv_errors dict, optimal_k)
    """
    admixture = get_admixture_executable()
    if admixture is None:
        return {}, 0

    ensure_dir(output_dir)

    # Copy input files to output directory (ADMIXTURE outputs to current directory)
    bed_file = input_prefix.parent / f"{input_prefix.name}.bed"
    bim_file = input_prefix.parent / f"{input_prefix.name}.bim"
    fam_file = input_prefix.parent / f"{input_prefix.name}.fam"

    work_bed = output_dir / bed_file.name
    work_bim = output_dir / bim_file.name
    work_fam = output_dir / fam_file.name

    shutil.copy(bed_file, work_bed)
    shutil.copy(bim_file, work_bim)
    shutil.copy(fam_file, work_fam)

    cv_errors: Dict[int, float] = {}
    k_range = range(min_k, max_k + 1)

    logger.info(f"Running ADMIXTURE for K={min_k} to K={max_k}...")

    # Run sequentially with progress bar
    for k in tqdm(k_range, desc="ADMIXTURE"):
        k_val, cv_error, success = run_admixture_k(
            work_bed, k, cv, threads, admixture, output_dir
        )

        if success:
            cv_errors[k_val] = cv_error
            logger.info(f"K={k_val}: CV error = {cv_error:.6f}")

            # Move output files to proper names
            q_file = output_dir / f"{work_bed.stem}.{k}.Q"
            p_file = output_dir / f"{work_bed.stem}.{k}.P"

            if q_file.exists():
                q_dest = output_dir / f"K{k}.Q"
                shutil.move(q_file, q_dest)

            if p_file.exists():
                p_dest = output_dir / f"K{k}.P"
                shutil.move(p_file, p_dest)
        else:
            logger.warning(f"ADMIXTURE failed for K={k}")

    # Clean up copied input files
    work_bed.unlink(missing_ok=True)
    work_bim.unlink(missing_ok=True)
    work_fam.unlink(missing_ok=True)

    # Find optimal K
    if cv_errors:
        optimal_k = min(cv_errors.keys(), key=lambda k: cv_errors[k])
        logger.info(f"Optimal K = {optimal_k} (CV error = {cv_errors[optimal_k]:.6f})")
    else:
        optimal_k = min_k

    return cv_errors, optimal_k


def save_cv_errors(
    cv_errors: Dict[int, float],
    output_file: Path,
) -> None:
    """Save cross-validation errors to file."""
    with open(output_file, "w") as f:
        f.write("K\tCV_Error\n")
        for k in sorted(cv_errors.keys()):
            f.write(f"{k}\t{cv_errors[k]:.6f}\n")

    logger.info(f"CV errors saved to: {output_file}")


def create_cv_plot(
    cv_errors: Dict[int, float],
    output_file: Path,
    optimal_k: int,
) -> None:
    """
    Create cross-validation error plot.

    Args:
        cv_errors: Dictionary of K -> CV error
        output_file: Output file path
        optimal_k: Optimal K value to highlight
    """
    fig, ax = plt.subplots(figsize=(10, 6))

    k_values = sorted(cv_errors.keys())
    errors = [cv_errors[k] for k in k_values]

    ax.plot(k_values, errors, "o-", color="steelblue", linewidth=2, markersize=8)

    # Highlight optimal K
    ax.scatter([optimal_k], [cv_errors[optimal_k]], color="red", s=150, zorder=5)
    ax.annotate(
        f"Optimal K={optimal_k}",
        xy=(optimal_k, cv_errors[optimal_k]),
        xytext=(optimal_k + 0.5, cv_errors[optimal_k]),
        fontsize=10,
        color="red",
    )

    ax.set_xlabel("K (Number of Ancestral Populations)", fontsize=12)
    ax.set_ylabel("Cross-Validation Error", fontsize=12)
    ax.set_title("ADMIXTURE Cross-Validation", fontsize=14)
    ax.set_xticks(k_values)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_file, dpi=150, bbox_inches="tight")
    plt.close()

    logger.info(f"CV plot saved to: {output_file}")


def load_admixture_q(q_file: Path) -> pd.DataFrame:
    """Load ADMIXTURE Q file."""
    return pd.read_csv(q_file, sep=r"\s+", header=None)


def create_admixture_barplot(
    q_df: pd.DataFrame,
    k: int,
    sample_ids: List[str],
    output_file: Path,
    populations: Optional[List[str]] = None,
    sort_by_population: bool = True,
    figsize: Tuple[int, int] = (20, 6),
) -> None:
    """
    Create ADMIXTURE ancestry proportion bar plot.

    Args:
        q_df: DataFrame with ancestry proportions
        k: K value
        sample_ids: List of sample IDs
        output_file: Output file path
        populations: Optional list of population labels
        sort_by_population: Whether to sort samples by population
        figsize: Figure size
    """
    fig, ax = plt.subplots(figsize=figsize)

    n_samples = len(q_df)

    # Create DataFrame with sample IDs
    q_df = q_df.copy()
    q_df["Sample"] = sample_ids
    if populations:
        q_df["Population"] = populations

    # Sort by population if requested
    if sort_by_population and populations:
        q_df = q_df.sort_values("Population")
        sample_order = q_df.index.tolist()
    else:
        sample_order = list(range(n_samples))

    # Define colors for ancestry components
    component_colors = [
        "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
        "#FFFF33", "#A65628", "#F781BF", "#999999", "#66C2A5",
        "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F",
    ]

    # Plot stacked bars
    x = np.arange(n_samples)
    bottom = np.zeros(n_samples)

    for i in range(k):
        values = q_df.iloc[sample_order, i].values
        ax.bar(
            x, values, bottom=bottom,
            color=component_colors[i % len(component_colors)],
            width=1.0,
            label=f"Ancestry {i+1}",
        )
        bottom += values

    # Add population labels if available
    if populations:
        # Find population boundaries
        sorted_pops = q_df.iloc[sample_order]["Population"].values
        boundaries = [0]
        for i in range(1, len(sorted_pops)):
            if sorted_pops[i] != sorted_pops[i-1]:
                boundaries.append(i)
        boundaries.append(len(sorted_pops))

        # Add population labels
        for i in range(len(boundaries) - 1):
            mid = (boundaries[i] + boundaries[i+1]) / 2
            pop_label = sorted_pops[boundaries[i]]
            ax.text(
                mid, -0.05, pop_label,
                ha="center", va="top",
                fontsize=6, rotation=90,
            )

    ax.set_xlim(-0.5, n_samples - 0.5)
    ax.set_ylim(0, 1)
    ax.set_ylabel("Ancestry Proportion", fontsize=12)
    ax.set_title(f"ADMIXTURE K={k}", fontsize=14)
    ax.set_xticks([])

    # Legend
    ax.legend(
        bbox_to_anchor=(1.02, 1),
        loc="upper left",
        fontsize=8,
    )

    plt.tight_layout()
    plt.savefig(output_file, dpi=150, bbox_inches="tight")
    plt.close()

    logger.info(f"Bar plot saved to: {output_file}")


def create_summary_table(
    output_dir: Path,
    fam_file: Path,
    k_values: List[int],
    population_file: Optional[Path] = None,
) -> None:
    """
    Create summary table of ancestry proportions.

    Args:
        output_dir: Directory with Q files
        fam_file: FAM file for sample IDs
        k_values: List of K values to include
        population_file: Optional population labels file
    """
    # Load sample information
    fam_df = read_fam_file(fam_file)
    sample_ids = fam_df["IID"].tolist()

    # Load population labels if available
    if population_file and population_file.exists():
        pop_map = read_population_file(population_file)
        populations = [pop_map.get(sid, "Unknown") for sid in sample_ids]
    else:
        populations = fam_df["FID"].tolist()

    # Create summary for each K
    for k in k_values:
        q_file = output_dir / f"K{k}.Q"
        if not q_file.exists():
            continue

        q_df = load_admixture_q(q_file)

        # Create summary DataFrame
        summary = pd.DataFrame({
            "Sample": sample_ids,
            "Population": populations,
        })

        for i in range(k):
            summary[f"Ancestry_{i+1}"] = q_df.iloc[:, i]

        # Save summary
        summary_file = output_dir / f"K{k}_summary.csv"
        summary.to_csv(summary_file, index=False)

        # Also save population means
        pop_means = summary.groupby("Population").mean(numeric_only=True)
        pop_means_file = output_dir / f"K{k}_population_means.csv"
        pop_means.to_csv(pop_means_file)


def run_admixture_pipeline(
    input_prefix: Path,
    output_dir: Path,
    min_k: int = 2,
    max_k: int = 15,
    cv: int = 5,
    threads: int = 4,
    population_file: Optional[Path] = None,
) -> bool:
    """
    Run the complete ADMIXTURE pipeline.

    Args:
        input_prefix: Input PLINK file prefix
        output_dir: Output directory
        min_k: Minimum K value
        max_k: Maximum K value
        cv: Cross-validation folds
        threads: Number of threads
        population_file: Optional population labels file

    Returns:
        True if successful
    """
    ensure_dir(output_dir)

    # Run ADMIXTURE for all K values
    cv_errors, optimal_k = run_admixture_range(
        input_prefix, output_dir, min_k, max_k, cv, threads
    )

    if not cv_errors:
        logger.error("ADMIXTURE failed for all K values")
        return False

    # Save CV errors
    save_cv_errors(cv_errors, output_dir / "cv_errors.txt")

    # Create CV plot
    create_cv_plot(cv_errors, output_dir / "cv_plot.png", optimal_k)

    # Load sample information
    fam_file = input_prefix.parent / f"{input_prefix.name}.fam"
    fam_df = read_fam_file(fam_file)
    sample_ids = fam_df["IID"].tolist()

    # Load population labels if available
    if population_file and population_file.exists():
        pop_map = read_population_file(population_file)
        populations = [pop_map.get(sid, "Unknown") for sid in sample_ids]
    else:
        populations = fam_df["FID"].tolist()

    # Create bar plots for each K
    logger.info("Creating bar plots...")

    for k in tqdm(cv_errors.keys(), desc="Bar plots"):
        q_file = output_dir / f"K{k}.Q"
        if q_file.exists():
            q_df = load_admixture_q(q_file)
            create_admixture_barplot(
                q_df, k, sample_ids,
                output_dir / f"admixture_plot_K{k}.png",
                populations=populations,
            )

    # Create summary tables
    create_summary_table(
        output_dir, fam_file,
        list(cv_errors.keys()),
        population_file,
    )

    logger.info("ADMIXTURE pipeline completed successfully")
    return True


def main() -> int:
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="ADMIXTURE analysis for ancestry proportions.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Basic ADMIXTURE analysis
    %(prog)s --input merged --output-dir results/admixture/

    # With custom K range
    %(prog)s --input merged --output-dir results/admixture/ \\
             --min-k 2 --max-k 8 --threads 8

Output:
    - K{k}.Q - Ancestry proportions for each K
    - K{k}.P - Allele frequencies for each K
    - cv_errors.txt - Cross-validation errors
    - cv_plot.png - CV error plot
    - admixture_plot_K{k}.png - Bar plots for each K
    - K{k}_summary.csv - Summary tables
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
        default=str(RESULTS_DIR / "admixture"),
        help=f"Output directory (default: {RESULTS_DIR / 'admixture'})",
    )

    parser.add_argument(
        "--min-k",
        type=int,
        default=ADMIXTURE_DEFAULTS["min_k"],
        help=f"Minimum K value (default: {ADMIXTURE_DEFAULTS['min_k']})",
    )

    parser.add_argument(
        "--max-k",
        type=int,
        default=ADMIXTURE_DEFAULTS["max_k"],
        help=f"Maximum K value (default: {ADMIXTURE_DEFAULTS['max_k']})",
    )

    parser.add_argument(
        "--cv",
        type=int,
        default=ADMIXTURE_DEFAULTS["cv_folds"],
        help=f"Cross-validation folds (default: {ADMIXTURE_DEFAULTS['cv_folds']})",
    )

    parser.add_argument(
        "--threads", "-j",
        type=int,
        default=ADMIXTURE_DEFAULTS["threads"],
        help=f"Number of threads (default: {ADMIXTURE_DEFAULTS['threads']})",
    )

    parser.add_argument(
        "--population-file",
        type=str,
        help="File mapping sample IDs to populations (optional)",
    )

    parser.add_argument(
        "--log-file",
        type=str,
        help="Path to log file",
    )

    args = parser.parse_args()

    # Set up logging
    log_file = Path(args.log_file) if args.log_file else LOGS_DIR / "run_admixture.log"
    setup_logging(log_file)

    logger.info("=" * 60)
    logger.info("ADMIXTURE Analysis")
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
    logger.info(f"K range: {args.min_k} to {args.max_k}")
    logger.info(f"CV folds: {args.cv}")
    logger.info(f"Threads: {args.threads}")

    # Run ADMIXTURE pipeline
    success = run_admixture_pipeline(
        input_prefix=input_prefix,
        output_dir=output_dir,
        min_k=args.min_k,
        max_k=args.max_k,
        cv=args.cv,
        threads=args.threads,
        population_file=population_file,
    )

    if success:
        logger.info("=" * 60)
        logger.info("ADMIXTURE completed successfully!")
        logger.info(f"Results saved to: {output_dir}")
        logger.info("=" * 60)
        return 0
    else:
        logger.error("ADMIXTURE failed")
        return 1


if __name__ == "__main__":
    sys.exit(main())
