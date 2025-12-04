#!/usr/bin/env python3
"""
Quality control filtering for genetic data.

This script performs comprehensive quality control on PLINK format genetic data,
including SNP and individual filtering, LD pruning, and relatedness filtering.

Features:
    - Filter SNPs by missingness (default: >5%)
    - Filter individuals by missingness (default: >10%)
    - Filter by minor allele frequency (MAF, default: >0.01)
    - Remove palindromic SNPs (A/T, G/C)
    - LD pruning (r² threshold, window size)
    - Identity-by-descent (IBD) filtering to remove related individuals
    - Generate QC report with statistics

Usage:
    python qc_filtering.py --input data/processed/sample \\
                           --output data/processed/sample_qc \\
                           --remove-palindromic

Example:
    # Basic QC with default parameters
    python qc_filtering.py --input sample --output sample_qc

    # Stringent QC with custom thresholds
    python qc_filtering.py --input sample --output sample_qc \\
                           --snp-missingness 0.02 --maf 0.05 --ld-r2 0.1
"""

import argparse
import logging
import subprocess
import sys
import tempfile
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from tqdm import tqdm

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from utils.config import (
    LOGS_DIR,
    QC_DEFAULTS,
    LD_DEFAULTS,
    RELATEDNESS_DEFAULTS,
    get_tool_path,
)
from utils.file_utils import (
    check_plink_files,
    ensure_dir,
    read_bim_file,
    count_plink_snps,
    count_plink_samples,
)
from utils.genetics_utils import is_palindromic


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


def run_plink_command(
    cmd: List[str],
    description: str = "Running PLINK",
) -> Tuple[bool, str, str]:
    """
    Run a PLINK command and capture output.

    Args:
        cmd: Command and arguments as list
        description: Description for logging

    Returns:
        Tuple of (success, stdout, stderr)
    """
    logger.info(f"{description}...")
    logger.debug(f"Command: {' '.join(cmd)}")

    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=False,
        )

        if result.returncode != 0:
            logger.error(f"PLINK error: {result.stderr}")
            return False, result.stdout, result.stderr

        return True, result.stdout, result.stderr

    except FileNotFoundError:
        logger.error(f"PLINK executable not found: {cmd[0]}")
        return False, "", "PLINK not found"
    except Exception as e:
        logger.error(f"Error running PLINK: {e}")
        return False, "", str(e)


def get_plink_executable() -> Optional[Path]:
    """Get PLINK executable path."""
    plink_path = get_tool_path("plink")
    if plink_path is None:
        logger.error("PLINK not found. Please install PLINK and ensure it's in PATH.")
        return None
    return plink_path


def identify_palindromic_snps(bim_file: Path) -> List[str]:
    """
    Identify palindromic SNPs (A/T, G/C) in a BIM file.

    Args:
        bim_file: Path to PLINK .bim file

    Returns:
        List of palindromic SNP IDs
    """
    logger.info("Identifying palindromic SNPs...")

    palindromic = []

    with open(bim_file, "r") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 6:
                snp_id = parts[1]
                a1 = parts[4].upper()
                a2 = parts[5].upper()

                if is_palindromic(a1, a2):
                    palindromic.append(snp_id)

    logger.info(f"Found {len(palindromic)} palindromic SNPs")
    return palindromic


def filter_snp_missingness(
    input_prefix: Path,
    output_prefix: Path,
    plink: Path,
    threshold: float = 0.05,
) -> Tuple[bool, int]:
    """
    Filter SNPs by missingness rate.

    Args:
        input_prefix: Input PLINK file prefix
        output_prefix: Output PLINK file prefix
        plink: Path to PLINK executable
        threshold: Maximum allowed missingness rate

    Returns:
        Tuple of (success, number of SNPs removed)
    """
    snps_before = count_plink_snps(input_prefix.parent / f"{input_prefix.name}.bim")

    cmd = [
        str(plink),
        "--bfile", str(input_prefix),
        "--geno", str(threshold),
        "--make-bed",
        "--out", str(output_prefix),
        "--allow-no-sex",
    ]

    success, _, _ = run_plink_command(cmd, f"Filtering SNPs with missingness > {threshold}")

    if success and check_plink_files(output_prefix):
        snps_after = count_plink_snps(output_prefix.parent / f"{output_prefix.name}.bim")
        removed = snps_before - snps_after
        logger.info(f"Removed {removed} SNPs due to high missingness")
        return True, removed

    return False, 0


def filter_ind_missingness(
    input_prefix: Path,
    output_prefix: Path,
    plink: Path,
    threshold: float = 0.10,
) -> Tuple[bool, int]:
    """
    Filter individuals by missingness rate.

    Args:
        input_prefix: Input PLINK file prefix
        output_prefix: Output PLINK file prefix
        plink: Path to PLINK executable
        threshold: Maximum allowed missingness rate

    Returns:
        Tuple of (success, number of individuals removed)
    """
    inds_before = count_plink_samples(input_prefix.parent / f"{input_prefix.name}.fam")

    cmd = [
        str(plink),
        "--bfile", str(input_prefix),
        "--mind", str(threshold),
        "--make-bed",
        "--out", str(output_prefix),
        "--allow-no-sex",
    ]

    success, _, _ = run_plink_command(cmd, f"Filtering individuals with missingness > {threshold}")

    if success and check_plink_files(output_prefix):
        inds_after = count_plink_samples(output_prefix.parent / f"{output_prefix.name}.fam")
        removed = inds_before - inds_after
        logger.info(f"Removed {removed} individuals due to high missingness")
        return True, removed

    return False, 0


def filter_maf(
    input_prefix: Path,
    output_prefix: Path,
    plink: Path,
    threshold: float = 0.01,
) -> Tuple[bool, int]:
    """
    Filter SNPs by minor allele frequency.

    Args:
        input_prefix: Input PLINK file prefix
        output_prefix: Output PLINK file prefix
        plink: Path to PLINK executable
        threshold: Minimum MAF threshold

    Returns:
        Tuple of (success, number of SNPs removed)
    """
    snps_before = count_plink_snps(input_prefix.parent / f"{input_prefix.name}.bim")

    cmd = [
        str(plink),
        "--bfile", str(input_prefix),
        "--maf", str(threshold),
        "--make-bed",
        "--out", str(output_prefix),
        "--allow-no-sex",
    ]

    success, _, _ = run_plink_command(cmd, f"Filtering SNPs with MAF < {threshold}")

    if success and check_plink_files(output_prefix):
        snps_after = count_plink_snps(output_prefix.parent / f"{output_prefix.name}.bim")
        removed = snps_before - snps_after
        logger.info(f"Removed {removed} SNPs due to low MAF")
        return True, removed

    return False, 0


def remove_palindromic_snps(
    input_prefix: Path,
    output_prefix: Path,
    plink: Path,
) -> Tuple[bool, int]:
    """
    Remove palindromic SNPs (A/T, G/C).

    Args:
        input_prefix: Input PLINK file prefix
        output_prefix: Output PLINK file prefix
        plink: Path to PLINK executable

    Returns:
        Tuple of (success, number of SNPs removed)
    """
    bim_file = input_prefix.parent / f"{input_prefix.name}.bim"
    palindromic = identify_palindromic_snps(bim_file)

    if not palindromic:
        # No palindromic SNPs, just copy files
        logger.info("No palindromic SNPs to remove")
        cmd = [
            str(plink),
            "--bfile", str(input_prefix),
            "--make-bed",
            "--out", str(output_prefix),
            "--allow-no-sex",
        ]
        success, _, _ = run_plink_command(cmd, "Copying files (no palindromic SNPs)")
        return success, 0

    # Write list of SNPs to exclude
    with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False) as f:
        exclude_file = Path(f.name)
        for snp_id in palindromic:
            f.write(f"{snp_id}\n")

    try:
        cmd = [
            str(plink),
            "--bfile", str(input_prefix),
            "--exclude", str(exclude_file),
            "--make-bed",
            "--out", str(output_prefix),
            "--allow-no-sex",
        ]

        success, _, _ = run_plink_command(cmd, "Removing palindromic SNPs")

        if success:
            logger.info(f"Removed {len(palindromic)} palindromic SNPs")
            return True, len(palindromic)

        return False, 0

    finally:
        exclude_file.unlink(missing_ok=True)


def ld_prune(
    input_prefix: Path,
    output_prefix: Path,
    plink: Path,
    window_kb: int = 50,
    step: int = 5,
    r2_threshold: float = 0.2,
) -> Tuple[bool, int]:
    """
    Perform LD pruning.

    Args:
        input_prefix: Input PLINK file prefix
        output_prefix: Output PLINK file prefix
        plink: Path to PLINK executable
        window_kb: Window size in kb
        step: Step size in SNPs
        r2_threshold: r² threshold

    Returns:
        Tuple of (success, number of SNPs removed)
    """
    snps_before = count_plink_snps(input_prefix.parent / f"{input_prefix.name}.bim")

    # Step 1: Generate list of SNPs to prune
    prune_prefix = output_prefix.parent / f"{output_prefix.name}_ldprune"

    cmd = [
        str(plink),
        "--bfile", str(input_prefix),
        "--indep-pairwise", str(window_kb), str(step), str(r2_threshold),
        "--out", str(prune_prefix),
        "--allow-no-sex",
    ]

    success, _, _ = run_plink_command(cmd, f"LD pruning (window={window_kb}kb, r²={r2_threshold})")

    if not success:
        return False, 0

    # Step 2: Extract pruned SNPs
    prune_in_file = prune_prefix.parent / f"{prune_prefix.name}.prune.in"

    if not prune_in_file.exists():
        logger.error("LD pruning output file not found")
        return False, 0

    cmd = [
        str(plink),
        "--bfile", str(input_prefix),
        "--extract", str(prune_in_file),
        "--make-bed",
        "--out", str(output_prefix),
        "--allow-no-sex",
    ]

    success, _, _ = run_plink_command(cmd, "Extracting LD-pruned SNPs")

    # Clean up intermediate files
    for suffix in [".prune.in", ".prune.out", ".log", ".nosex"]:
        temp_file = prune_prefix.parent / f"{prune_prefix.name}{suffix}"
        temp_file.unlink(missing_ok=True)

    if success and check_plink_files(output_prefix):
        snps_after = count_plink_snps(output_prefix.parent / f"{output_prefix.name}.bim")
        removed = snps_before - snps_after
        logger.info(f"Removed {removed} SNPs through LD pruning")
        return True, removed

    return False, 0


def filter_relatedness(
    input_prefix: Path,
    output_prefix: Path,
    plink: Path,
    pi_hat_threshold: float = 0.2,
) -> Tuple[bool, int]:
    """
    Filter related individuals using IBD (PI_HAT).

    Args:
        input_prefix: Input PLINK file prefix
        output_prefix: Output PLINK file prefix
        plink: Path to PLINK executable
        pi_hat_threshold: PI_HAT threshold for relatedness

    Returns:
        Tuple of (success, number of individuals removed)
    """
    inds_before = count_plink_samples(input_prefix.parent / f"{input_prefix.name}.fam")

    # Calculate IBD
    ibd_prefix = output_prefix.parent / f"{output_prefix.name}_ibd"

    cmd = [
        str(plink),
        "--bfile", str(input_prefix),
        "--genome",
        "--out", str(ibd_prefix),
        "--allow-no-sex",
    ]

    success, _, _ = run_plink_command(cmd, "Calculating IBD (genome-wide identity)")

    if not success:
        return False, 0

    # Parse genome file and identify related individuals
    genome_file = ibd_prefix.parent / f"{ibd_prefix.name}.genome"

    if not genome_file.exists():
        logger.warning("No genome file generated - skipping relatedness filter")
        # Just copy the files
        cmd = [
            str(plink),
            "--bfile", str(input_prefix),
            "--make-bed",
            "--out", str(output_prefix),
            "--allow-no-sex",
        ]
        success, _, _ = run_plink_command(cmd, "Copying files (no related individuals)")
        return success, 0

    # Find pairs with PI_HAT > threshold
    related_pairs = []
    with open(genome_file, "r") as f:
        header = f.readline()  # Skip header
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 10:
                pi_hat = float(parts[9])  # PI_HAT column
                if pi_hat > pi_hat_threshold:
                    fid1, iid1 = parts[0], parts[1]
                    fid2, iid2 = parts[2], parts[3]
                    related_pairs.append(((fid1, iid1), (fid2, iid2), pi_hat))

    if not related_pairs:
        logger.info("No related individuals found")
        cmd = [
            str(plink),
            "--bfile", str(input_prefix),
            "--make-bed",
            "--out", str(output_prefix),
            "--allow-no-sex",
        ]
        success, _, _ = run_plink_command(cmd, "Copying files (no related individuals)")

        # Clean up
        for suffix in [".genome", ".log", ".nosex"]:
            temp_file = ibd_prefix.parent / f"{ibd_prefix.name}{suffix}"
            temp_file.unlink(missing_ok=True)

        return success, 0

    # Determine which individuals to remove (greedy approach: remove one from each pair)
    # In case of multiple relationships, remove the individual involved in most pairs
    individual_counts: Dict[Tuple[str, str], int] = {}
    for (ind1, ind2, _) in related_pairs:
        individual_counts[ind1] = individual_counts.get(ind1, 0) + 1
        individual_counts[ind2] = individual_counts.get(ind2, 0) + 1

    to_remove = set()
    for (ind1, ind2, _) in related_pairs:
        if ind1 in to_remove or ind2 in to_remove:
            continue
        # Remove the one with more relationships
        if individual_counts[ind1] >= individual_counts[ind2]:
            to_remove.add(ind1)
        else:
            to_remove.add(ind2)

    if to_remove:
        # Write list of individuals to remove
        with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False) as f:
            remove_file = Path(f.name)
            for fid, iid in to_remove:
                f.write(f"{fid}\t{iid}\n")

        try:
            cmd = [
                str(plink),
                "--bfile", str(input_prefix),
                "--remove", str(remove_file),
                "--make-bed",
                "--out", str(output_prefix),
                "--allow-no-sex",
            ]

            success, _, _ = run_plink_command(cmd, "Removing related individuals")

            if success:
                logger.info(f"Removed {len(to_remove)} related individuals (PI_HAT > {pi_hat_threshold})")

        finally:
            remove_file.unlink(missing_ok=True)
    else:
        cmd = [
            str(plink),
            "--bfile", str(input_prefix),
            "--make-bed",
            "--out", str(output_prefix),
            "--allow-no-sex",
        ]
        success, _, _ = run_plink_command(cmd, "Copying files")

    # Clean up
    for suffix in [".genome", ".log", ".nosex"]:
        temp_file = ibd_prefix.parent / f"{ibd_prefix.name}{suffix}"
        temp_file.unlink(missing_ok=True)

    if success and check_plink_files(output_prefix):
        inds_after = count_plink_samples(output_prefix.parent / f"{output_prefix.name}.fam")
        removed = inds_before - inds_after
        return True, removed

    return False, 0


def write_qc_report(
    stats: Dict,
    output_path: Path,
) -> None:
    """
    Write QC statistics report.

    Args:
        stats: Dictionary of QC statistics
        output_path: Output file path
    """
    with open(output_path, "w") as f:
        f.write("=" * 60 + "\n")
        f.write("Quality Control Report\n")
        f.write("=" * 60 + "\n")
        f.write(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write("\n")

        f.write("Input Data:\n")
        f.write(f"  Input file: {stats.get('input_prefix', 'N/A')}\n")
        f.write(f"  Initial SNPs: {stats.get('initial_snps', 'N/A')}\n")
        f.write(f"  Initial samples: {stats.get('initial_samples', 'N/A')}\n")
        f.write("\n")

        f.write("Parameters:\n")
        f.write(f"  SNP missingness threshold: {stats.get('snp_missingness_threshold', 'N/A')}\n")
        f.write(f"  Individual missingness threshold: {stats.get('ind_missingness_threshold', 'N/A')}\n")
        f.write(f"  MAF threshold: {stats.get('maf_threshold', 'N/A')}\n")
        f.write(f"  LD window (kb): {stats.get('ld_window', 'N/A')}\n")
        f.write(f"  LD r² threshold: {stats.get('ld_r2', 'N/A')}\n")
        f.write(f"  PI_HAT threshold: {stats.get('pi_hat', 'N/A')}\n")
        f.write(f"  Remove palindromic: {stats.get('remove_palindromic', False)}\n")
        f.write("\n")

        f.write("QC Steps:\n")
        f.write(f"  SNPs removed (missingness): {stats.get('snps_removed_missingness', 0)}\n")
        f.write(f"  Individuals removed (missingness): {stats.get('inds_removed_missingness', 0)}\n")
        f.write(f"  SNPs removed (MAF): {stats.get('snps_removed_maf', 0)}\n")
        f.write(f"  SNPs removed (palindromic): {stats.get('snps_removed_palindromic', 0)}\n")
        f.write(f"  SNPs removed (LD pruning): {stats.get('snps_removed_ld', 0)}\n")
        f.write(f"  Individuals removed (relatedness): {stats.get('inds_removed_relatedness', 0)}\n")
        f.write("\n")

        f.write("Output Data:\n")
        f.write(f"  Output file: {stats.get('output_prefix', 'N/A')}\n")
        f.write(f"  Final SNPs: {stats.get('final_snps', 'N/A')}\n")
        f.write(f"  Final samples: {stats.get('final_samples', 'N/A')}\n")
        f.write("\n")

        total_snps_removed = (
            stats.get('snps_removed_missingness', 0) +
            stats.get('snps_removed_maf', 0) +
            stats.get('snps_removed_palindromic', 0) +
            stats.get('snps_removed_ld', 0)
        )
        total_inds_removed = (
            stats.get('inds_removed_missingness', 0) +
            stats.get('inds_removed_relatedness', 0)
        )

        f.write("Summary:\n")
        f.write(f"  Total SNPs removed: {total_snps_removed}\n")
        f.write(f"  Total individuals removed: {total_inds_removed}\n")
        f.write(f"  SNP retention rate: {stats.get('final_snps', 0) / max(stats.get('initial_snps', 1), 1) * 100:.1f}%\n")
        f.write(f"  Sample retention rate: {stats.get('final_samples', 0) / max(stats.get('initial_samples', 1), 1) * 100:.1f}%\n")
        f.write("=" * 60 + "\n")

    logger.info(f"QC report written to: {output_path}")


def run_qc_pipeline(
    input_prefix: Path,
    output_prefix: Path,
    snp_missingness: float = 0.05,
    ind_missingness: float = 0.10,
    maf: float = 0.01,
    ld_window: int = 50,
    ld_r2: float = 0.2,
    pi_hat: float = 0.2,
    remove_palindromic: bool = False,
) -> Tuple[bool, Dict]:
    """
    Run the complete QC pipeline.

    Args:
        input_prefix: Input PLINK file prefix
        output_prefix: Output PLINK file prefix
        snp_missingness: Maximum SNP missingness rate
        ind_missingness: Maximum individual missingness rate
        maf: Minimum MAF threshold
        ld_window: LD pruning window size (kb)
        ld_r2: LD r² threshold
        pi_hat: PI_HAT threshold for relatedness
        remove_palindromic: Whether to remove palindromic SNPs

    Returns:
        Tuple of (success, statistics dict)
    """
    plink = get_plink_executable()
    if plink is None:
        return False, {}

    # Ensure output directory exists
    ensure_dir(output_prefix.parent)

    # Initialize statistics
    stats = {
        "input_prefix": str(input_prefix),
        "output_prefix": str(output_prefix),
        "initial_snps": count_plink_snps(input_prefix.parent / f"{input_prefix.name}.bim"),
        "initial_samples": count_plink_samples(input_prefix.parent / f"{input_prefix.name}.fam"),
        "snp_missingness_threshold": snp_missingness,
        "ind_missingness_threshold": ind_missingness,
        "maf_threshold": maf,
        "ld_window": ld_window,
        "ld_r2": ld_r2,
        "pi_hat": pi_hat,
        "remove_palindromic": remove_palindromic,
    }

    logger.info(f"Starting QC pipeline")
    logger.info(f"  Initial SNPs: {stats['initial_snps']}")
    logger.info(f"  Initial samples: {stats['initial_samples']}")

    # Use temporary directory for intermediate files
    temp_dir = output_prefix.parent / "temp_qc"
    ensure_dir(temp_dir)

    current_prefix = input_prefix
    step = 0

    try:
        # Step 1: Filter SNP missingness
        step += 1
        next_prefix = temp_dir / f"step{step}_snp_miss"
        success, removed = filter_snp_missingness(current_prefix, next_prefix, plink, snp_missingness)
        if not success:
            return False, stats
        stats["snps_removed_missingness"] = removed
        current_prefix = next_prefix

        # Step 2: Filter individual missingness
        step += 1
        next_prefix = temp_dir / f"step{step}_ind_miss"
        success, removed = filter_ind_missingness(current_prefix, next_prefix, plink, ind_missingness)
        if not success:
            return False, stats
        stats["inds_removed_missingness"] = removed
        current_prefix = next_prefix

        # Step 3: Filter MAF
        step += 1
        next_prefix = temp_dir / f"step{step}_maf"
        success, removed = filter_maf(current_prefix, next_prefix, plink, maf)
        if not success:
            return False, stats
        stats["snps_removed_maf"] = removed
        current_prefix = next_prefix

        # Step 4: Remove palindromic SNPs (optional)
        if remove_palindromic:
            step += 1
            next_prefix = temp_dir / f"step{step}_palindromic"
            success, removed = remove_palindromic_snps(current_prefix, next_prefix, plink)
            if not success:
                return False, stats
            stats["snps_removed_palindromic"] = removed
            current_prefix = next_prefix
        else:
            stats["snps_removed_palindromic"] = 0

        # Step 5: LD pruning
        step += 1
        next_prefix = temp_dir / f"step{step}_ld"
        success, removed = ld_prune(current_prefix, next_prefix, plink, ld_window, 5, ld_r2)
        if not success:
            return False, stats
        stats["snps_removed_ld"] = removed
        current_prefix = next_prefix

        # Step 6: Relatedness filtering
        step += 1
        success, removed = filter_relatedness(current_prefix, output_prefix, plink, pi_hat)
        if not success:
            return False, stats
        stats["inds_removed_relatedness"] = removed

        # Final statistics
        stats["final_snps"] = count_plink_snps(output_prefix.parent / f"{output_prefix.name}.bim")
        stats["final_samples"] = count_plink_samples(output_prefix.parent / f"{output_prefix.name}.fam")

        logger.info(f"QC pipeline completed")
        logger.info(f"  Final SNPs: {stats['final_snps']}")
        logger.info(f"  Final samples: {stats['final_samples']}")

        return True, stats

    finally:
        # Clean up temporary files
        import shutil
        if temp_dir.exists():
            shutil.rmtree(temp_dir)

        # Clean up PLINK log files
        for log_file in output_prefix.parent.glob("*.log"):
            log_file.unlink(missing_ok=True)
        for nosex_file in output_prefix.parent.glob("*.nosex"):
            nosex_file.unlink(missing_ok=True)


def main() -> int:
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Quality control filtering for genetic data.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Basic QC with default parameters
    %(prog)s --input sample --output sample_qc

    # Stringent QC
    %(prog)s --input sample --output sample_qc \\
             --snp-missingness 0.02 --maf 0.05 --ld-r2 0.1

    # Include palindromic SNP removal
    %(prog)s --input sample --output sample_qc --remove-palindromic

Output:
    - {output}.bed/bim/fam - Filtered PLINK files
    - {output}_qc_report.txt - QC statistics report
        """,
    )

    parser.add_argument(
        "--input", "-i",
        type=str,
        required=True,
        help="Input PLINK file prefix (required)",
    )

    parser.add_argument(
        "--output", "-o",
        type=str,
        required=True,
        help="Output PLINK file prefix (required)",
    )

    parser.add_argument(
        "--snp-missingness",
        type=float,
        default=QC_DEFAULTS["snp_missingness"],
        help=f"Maximum SNP missingness rate (default: {QC_DEFAULTS['snp_missingness']})",
    )

    parser.add_argument(
        "--ind-missingness",
        type=float,
        default=QC_DEFAULTS["ind_missingness"],
        help=f"Maximum individual missingness rate (default: {QC_DEFAULTS['ind_missingness']})",
    )

    parser.add_argument(
        "--maf",
        type=float,
        default=QC_DEFAULTS["maf"],
        help=f"Minimum minor allele frequency (default: {QC_DEFAULTS['maf']})",
    )

    parser.add_argument(
        "--ld-window",
        type=int,
        default=int(LD_DEFAULTS["window_kb"]),
        help=f"LD pruning window in kb (default: {int(LD_DEFAULTS['window_kb'])})",
    )

    parser.add_argument(
        "--ld-r2",
        type=float,
        default=LD_DEFAULTS["r2_threshold"],
        help=f"LD r² threshold (default: {LD_DEFAULTS['r2_threshold']})",
    )

    parser.add_argument(
        "--pi-hat",
        type=float,
        default=RELATEDNESS_DEFAULTS["pi_hat"],
        help=f"PI_HAT threshold for relatedness (default: {RELATEDNESS_DEFAULTS['pi_hat']})",
    )

    parser.add_argument(
        "--remove-palindromic",
        action="store_true",
        help="Remove palindromic SNPs (A/T, G/C)",
    )

    parser.add_argument(
        "--log-file",
        type=str,
        help="Path to log file",
    )

    args = parser.parse_args()

    # Set up logging
    log_file = Path(args.log_file) if args.log_file else LOGS_DIR / "qc_filtering.log"
    setup_logging(log_file)

    logger.info("=" * 60)
    logger.info("Quality Control Filtering")
    logger.info("=" * 60)

    # Validate input
    input_prefix = Path(args.input).resolve()
    if not check_plink_files(input_prefix):
        logger.error(f"Input PLINK files not found: {input_prefix}")
        return 1

    output_prefix = Path(args.output).resolve()

    # Run QC pipeline
    success, stats = run_qc_pipeline(
        input_prefix=input_prefix,
        output_prefix=output_prefix,
        snp_missingness=args.snp_missingness,
        ind_missingness=args.ind_missingness,
        maf=args.maf,
        ld_window=args.ld_window,
        ld_r2=args.ld_r2,
        pi_hat=args.pi_hat,
        remove_palindromic=args.remove_palindromic,
    )

    if success:
        # Write QC report
        report_path = output_prefix.parent / f"{output_prefix.name}_qc_report.txt"
        write_qc_report(stats, report_path)

        logger.info("=" * 60)
        logger.info("QC completed successfully!")
        logger.info("=" * 60)
        return 0
    else:
        logger.error("QC pipeline failed")
        return 1


if __name__ == "__main__":
    sys.exit(main())
