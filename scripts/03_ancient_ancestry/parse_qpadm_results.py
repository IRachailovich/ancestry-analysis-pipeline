#!/usr/bin/env python3
"""
Parse qpAdm output files and generate summary reports.

This script parses qpAdm output files to extract ancestry coefficients,
standard errors, p-values, and generates visualizations and summary tables.

Features:
    - Extract ancestry coefficients from qpAdm output files
    - Calculate standard errors
    - Extract p-values for model fit
    - Generate summary tables
    - Create visualization of ancient ancestry proportions
    - Compare multiple models

Usage:
    python parse_qpadm_results.py --input-dir results/qpadm/ \\
                                  --output-dir results/qpadm/

Example:
    python parse_qpadm_results.py --input-dir results/qpadm/
"""

import argparse
import logging
import re
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from utils.config import LOGS_DIR, RESULTS_DIR
from utils.file_utils import ensure_dir


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


def parse_qpadm_output(output_file: Path) -> Optional[Dict]:
    """
    Parse a single qpAdm output file.

    Args:
        output_file: Path to qpAdm output file

    Returns:
        Dictionary with parsed results, or None if parsing failed
    """
    if not output_file.exists():
        logger.warning(f"Output file not found: {output_file}")
        return None

    results = {
        "model": output_file.stem.replace("_output", ""),
        "file": str(output_file),
        "sources": [],
        "coefficients": [],
        "std_errors": [],
        "p_value": None,
        "chisq": None,
        "df": None,
        "target": None,
        "valid": False,
    }

    with open(output_file, "r") as f:
        content = f.read()

    lines = content.split("\n")

    # Parse target population
    target_match = re.search(r"left pops:\s*\n\s*(\S+)", content)
    if target_match:
        results["target"] = target_match.group(1)

    # Parse source populations
    # Look for lines like "source1: weight"
    in_coefficients = False
    coefficient_lines = []

    for i, line in enumerate(lines):
        line = line.strip()

        # Look for best coefficients section
        if "best coefficients" in line.lower():
            in_coefficients = True
            continue

        if in_coefficients:
            # Parse coefficient lines
            # Format varies, but typically: population weight std_err
            parts = line.split()

            if len(parts) >= 2:
                try:
                    # Try to parse as: pop weight [std_err]
                    pop = parts[0]
                    weight = float(parts[1])

                    # Skip if this looks like a header or non-data line
                    if pop.lower() in ["fixed", "tail", "chisq", "df", "pvalue"]:
                        in_coefficients = False
                        continue

                    results["sources"].append(pop)
                    results["coefficients"].append(weight)

                    if len(parts) >= 3:
                        try:
                            std_err = float(parts[2])
                            results["std_errors"].append(std_err)
                        except ValueError:
                            results["std_errors"].append(None)
                    else:
                        results["std_errors"].append(None)

                except ValueError:
                    # End of coefficients section
                    in_coefficients = False

            elif not parts:
                # Empty line might end the section
                in_coefficients = False

    # Parse p-value (tail probability)
    pvalue_patterns = [
        r"tail\s*prob[^:]*:\s*([\d.e+-]+)",
        r"pvalue[^:]*:\s*([\d.e+-]+)",
        r"p[- ]?value[^:]*:\s*([\d.e+-]+)",
    ]

    for pattern in pvalue_patterns:
        match = re.search(pattern, content, re.IGNORECASE)
        if match:
            try:
                results["p_value"] = float(match.group(1))
                break
            except ValueError:
                pass

    # Parse chi-squared
    chisq_match = re.search(r"chisq[^:]*:\s*([\d.e+-]+)", content, re.IGNORECASE)
    if chisq_match:
        try:
            results["chisq"] = float(chisq_match.group(1))
        except ValueError:
            pass

    # Parse degrees of freedom
    df_match = re.search(r"\bdf[^:]*:\s*(\d+)", content, re.IGNORECASE)
    if df_match:
        try:
            results["df"] = int(df_match.group(1))
        except ValueError:
            pass

    # Check if results are valid
    if results["coefficients"]:
        results["valid"] = True

        # Fill in missing standard errors
        while len(results["std_errors"]) < len(results["coefficients"]):
            results["std_errors"].append(None)

    # Alternative parsing for different qpAdm output formats
    if not results["valid"]:
        # Try parsing format: "pop: coeff stderr"
        pop_pattern = r"(\S+):\s+([\d.e+-]+)\s+([\d.e+-]+)"
        matches = re.findall(pop_pattern, content)

        if matches:
            for pop, coeff, stderr in matches:
                if pop.lower() not in ["fixed", "tail", "chisq", "df", "pvalue"]:
                    results["sources"].append(pop)
                    try:
                        results["coefficients"].append(float(coeff))
                        results["std_errors"].append(float(stderr))
                    except ValueError:
                        pass

            if results["coefficients"]:
                results["valid"] = True

    return results


def parse_all_outputs(input_dir: Path) -> List[Dict]:
    """
    Parse all qpAdm output files in a directory.

    Args:
        input_dir: Directory containing qpAdm output files

    Returns:
        List of parsed result dictionaries
    """
    results = []

    # Find all output files
    output_files = list(input_dir.glob("*_output.txt"))

    if not output_files:
        logger.warning(f"No output files found in {input_dir}")
        return results

    logger.info(f"Found {len(output_files)} output file(s)")

    for output_file in sorted(output_files):
        logger.info(f"Parsing: {output_file.name}")
        parsed = parse_qpadm_output(output_file)

        if parsed and parsed["valid"]:
            results.append(parsed)
            logger.info(f"  Sources: {', '.join(parsed['sources'])}")
            logger.info(f"  Coefficients: {parsed['coefficients']}")
            if parsed["p_value"]:
                logger.info(f"  P-value: {parsed['p_value']:.4g}")
        else:
            logger.warning(f"  Failed to parse or invalid results")

    return results


def create_summary_table(
    results: List[Dict],
    output_file: Path,
) -> pd.DataFrame:
    """
    Create summary table from parsed results.

    Args:
        results: List of parsed result dictionaries
        output_file: Output CSV file path

    Returns:
        Summary DataFrame
    """
    rows = []

    for result in results:
        if not result["valid"]:
            continue

        row = {
            "Model": result["model"],
            "Target": result["target"],
            "P_value": result["p_value"],
            "Chi_sq": result["chisq"],
            "DF": result["df"],
        }

        # Add source columns
        for i, (source, coeff, stderr) in enumerate(
            zip(result["sources"], result["coefficients"], result["std_errors"])
        ):
            row[f"Source_{i+1}"] = source
            row[f"Coeff_{i+1}"] = coeff
            row[f"SE_{i+1}"] = stderr

        rows.append(row)

    df = pd.DataFrame(rows)

    if not df.empty:
        df.to_csv(output_file, index=False)
        logger.info(f"Summary table saved to: {output_file}")

    return df


def create_coefficient_plot(
    results: List[Dict],
    output_file: Path,
) -> None:
    """
    Create bar plot of ancestry coefficients.

    Args:
        results: List of parsed result dictionaries
        output_file: Output image file path
    """
    if not results:
        logger.warning("No valid results to plot")
        return

    # Filter valid results
    valid_results = [r for r in results if r["valid"] and r["coefficients"]]

    if not valid_results:
        logger.warning("No valid results with coefficients")
        return

    n_models = len(valid_results)
    fig, axes = plt.subplots(1, n_models, figsize=(6 * n_models, 6), squeeze=False)

    colors = ["#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00"]

    for idx, result in enumerate(valid_results):
        ax = axes[0, idx]

        sources = result["sources"]
        coefficients = result["coefficients"]
        std_errors = result["std_errors"]

        # Replace None with 0 for plotting
        std_errors = [se if se is not None else 0 for se in std_errors]

        x = np.arange(len(sources))
        bar_colors = [colors[i % len(colors)] for i in range(len(sources))]

        bars = ax.bar(x, coefficients, color=bar_colors, yerr=std_errors, capsize=5)

        ax.set_xticks(x)
        ax.set_xticklabels(sources, rotation=45, ha="right")
        ax.set_ylabel("Ancestry Proportion")
        ax.set_title(f"Model: {result['model']}")
        ax.set_ylim(0, 1)

        # Add value labels on bars
        for bar, coeff in zip(bars, coefficients):
            ax.text(
                bar.get_x() + bar.get_width() / 2,
                bar.get_height() + 0.02,
                f"{coeff:.2f}",
                ha="center",
                va="bottom",
                fontsize=9,
            )

        # Add p-value annotation
        if result["p_value"] is not None:
            pval_text = f"p = {result['p_value']:.4g}"
            ax.text(
                0.95, 0.95, pval_text,
                transform=ax.transAxes,
                ha="right", va="top",
                fontsize=10,
                bbox=dict(boxstyle="round", facecolor="white", alpha=0.8),
            )

    plt.tight_layout()
    plt.savefig(output_file, dpi=150, bbox_inches="tight")
    plt.close()

    logger.info(f"Coefficient plot saved to: {output_file}")


def create_model_comparison_plot(
    results: List[Dict],
    output_file: Path,
) -> None:
    """
    Create comparison plot across models.

    Args:
        results: List of parsed result dictionaries
        output_file: Output image file path
    """
    valid_results = [r for r in results if r["valid"] and r["p_value"] is not None]

    if len(valid_results) < 2:
        logger.info("Not enough models for comparison plot")
        return

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

    models = [r["model"] for r in valid_results]
    p_values = [r["p_value"] for r in valid_results]

    # P-value comparison
    colors = ["green" if p > 0.05 else "red" for p in p_values]
    ax1.barh(models, p_values, color=colors)
    ax1.axvline(x=0.05, color="black", linestyle="--", label="p = 0.05")
    ax1.set_xlabel("P-value")
    ax1.set_title("Model Fit Comparison")
    ax1.legend()

    # Coefficient comparison (first coefficient only for simplicity)
    first_coeffs = [r["coefficients"][0] if r["coefficients"] else 0 for r in valid_results]
    first_sources = [r["sources"][0] if r["sources"] else "Unknown" for r in valid_results]

    ax2.barh(models, first_coeffs, color="steelblue")
    ax2.set_xlabel(f"Coefficient ({first_sources[0] if first_sources else 'Source 1'})")
    ax2.set_title("First Source Coefficient Comparison")

    plt.tight_layout()
    plt.savefig(output_file, dpi=150, bbox_inches="tight")
    plt.close()

    logger.info(f"Model comparison plot saved to: {output_file}")


def generate_text_report(
    results: List[Dict],
    output_file: Path,
) -> None:
    """
    Generate detailed text report.

    Args:
        results: List of parsed result dictionaries
        output_file: Output text file path
    """
    with open(output_file, "w") as f:
        f.write("=" * 70 + "\n")
        f.write("                    qpAdm Results Summary                    \n")
        f.write("=" * 70 + "\n\n")

        for result in results:
            f.write(f"Model: {result['model']}\n")
            f.write("-" * 50 + "\n")

            if result["target"]:
                f.write(f"Target: {result['target']}\n")

            if result["valid"]:
                f.write("\nAncestry Proportions:\n")
                for source, coeff, stderr in zip(
                    result["sources"], result["coefficients"], result["std_errors"]
                ):
                    se_str = f" ± {stderr:.4f}" if stderr else ""
                    f.write(f"  {source}: {coeff:.4f}{se_str}\n")

                f.write("\nModel Fit:\n")
                if result["p_value"] is not None:
                    f.write(f"  P-value: {result['p_value']:.4g}\n")
                    if result["p_value"] < 0.05:
                        f.write("  ⚠ Warning: P < 0.05 suggests poor model fit\n")
                    else:
                        f.write("  ✓ Good model fit (P ≥ 0.05)\n")

                if result["chisq"] is not None:
                    f.write(f"  Chi-squared: {result['chisq']:.4f}\n")
                if result["df"] is not None:
                    f.write(f"  Degrees of freedom: {result['df']}\n")
            else:
                f.write("  ⚠ Failed to parse results\n")

            f.write("\n" + "=" * 70 + "\n\n")

        # Summary
        valid_models = [r for r in results if r["valid"]]
        good_fit_models = [r for r in valid_models if r["p_value"] and r["p_value"] >= 0.05]

        f.write("SUMMARY\n")
        f.write("-" * 50 + "\n")
        f.write(f"Total models: {len(results)}\n")
        f.write(f"Valid models: {len(valid_models)}\n")
        f.write(f"Good fit models (p ≥ 0.05): {len(good_fit_models)}\n")

        if good_fit_models:
            f.write("\nRecommended models:\n")
            for model in good_fit_models:
                f.write(f"  - {model['model']} (p = {model['p_value']:.4g})\n")

    logger.info(f"Text report saved to: {output_file}")


def run_parser_pipeline(
    input_dir: Path,
    output_dir: Path,
) -> bool:
    """
    Run the complete qpAdm results parsing pipeline.

    Args:
        input_dir: Directory with qpAdm output files
        output_dir: Output directory for parsed results

    Returns:
        True if successful
    """
    ensure_dir(output_dir)

    # Parse all output files
    results = parse_all_outputs(input_dir)

    if not results:
        logger.error("No valid results parsed")
        return False

    # Create summary table
    create_summary_table(results, output_dir / "qpadm_summary.csv")

    # Create visualizations
    create_coefficient_plot(results, output_dir / "qpadm_plot.png")
    create_model_comparison_plot(results, output_dir / "qpadm_comparison.png")

    # Generate text report
    generate_text_report(results, output_dir / "qpadm_report.txt")

    logger.info("qpAdm results parsing completed")
    return True


def main() -> int:
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Parse qpAdm output files and generate summary reports.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Parse results in default directory
    %(prog)s --input-dir results/qpadm/

    # Specify custom output directory
    %(prog)s --input-dir results/qpadm/ --output-dir results/qpadm/parsed/

Output:
    - qpadm_summary.csv - Summary table of all models
    - qpadm_plot.png - Ancestry coefficient bar plots
    - qpadm_comparison.png - Model comparison plot
    - qpadm_report.txt - Detailed text report
        """,
    )

    parser.add_argument(
        "--input-dir", "-i",
        type=str,
        required=True,
        help="Directory with qpAdm output files (required)",
    )

    parser.add_argument(
        "--output-dir", "-o",
        type=str,
        help="Output directory (default: same as input-dir)",
    )

    parser.add_argument(
        "--log-file",
        type=str,
        help="Path to log file",
    )

    args = parser.parse_args()

    # Set up logging
    log_file = Path(args.log_file) if args.log_file else LOGS_DIR / "parse_qpadm_results.log"
    setup_logging(log_file)

    logger.info("=" * 60)
    logger.info("qpAdm Results Parser")
    logger.info("=" * 60)

    # Set directories
    input_dir = Path(args.input_dir).resolve()
    output_dir = Path(args.output_dir).resolve() if args.output_dir else input_dir

    if not input_dir.exists():
        logger.error(f"Input directory not found: {input_dir}")
        return 1

    logger.info(f"Input directory: {input_dir}")
    logger.info(f"Output directory: {output_dir}")

    # Run parser pipeline
    success = run_parser_pipeline(input_dir, output_dir)

    if success:
        logger.info("=" * 60)
        logger.info("qpAdm results parsing completed successfully!")
        logger.info(f"Results saved to: {output_dir}")
        logger.info("=" * 60)
        return 0
    else:
        logger.error("qpAdm results parsing failed")
        return 1


if __name__ == "__main__":
    sys.exit(main())
