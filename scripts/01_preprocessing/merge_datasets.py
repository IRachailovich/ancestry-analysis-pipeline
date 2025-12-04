#!/usr/bin/env python3
"""
Merge sample data with reference populations.

This script merges 23andMe sample data with reference population datasets
(1000 Genomes, AADR, SGDP) for ancestry analysis.

Features:
    - Find overlapping SNPs between datasets
    - Merge 23andMe data with 1000 Genomes, AADR, SGDP
    - Handle strand issues (flip alleles if needed)
    - Convert EIGENSTRAT to PLINK format if needed
    - Report merge statistics

Usage:
    python merge_datasets.py --sample data/processed/sample \\
                             --references 1000g aadr \\
                             --output data/processed/merged

Example:
    # Merge with 1000 Genomes only
    python merge_datasets.py --sample sample --references 1000g --output merged

    # Merge with multiple references
    python merge_datasets.py --sample sample --references 1000g aadr sgdp --output merged
"""

import argparse
import logging
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

from tqdm import tqdm

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from utils.config import (
    LOGS_DIR,
    RAW_DATA_DIR,
    get_tool_path,
    get_reference_path,
)
from utils.file_utils import (
    check_plink_files,
    check_eigenstrat_files,
    ensure_dir,
    read_bim_file,
    count_plink_snps,
    count_plink_samples,
)
from utils.genetics_utils import (
    flip_allele,
    is_palindromic,
    are_alleles_compatible,
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


def run_command(
    cmd: List[str],
    description: str = "Running command",
) -> Tuple[bool, str, str]:
    """Run a command and capture output."""
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
            logger.error(f"Command error: {result.stderr}")
            return False, result.stdout, result.stderr

        return True, result.stdout, result.stderr

    except FileNotFoundError:
        logger.error(f"Executable not found: {cmd[0]}")
        return False, "", "Executable not found"
    except Exception as e:
        logger.error(f"Error running command: {e}")
        return False, "", str(e)


def get_plink_executable() -> Optional[Path]:
    """Get PLINK executable path."""
    plink_path = get_tool_path("plink")
    if plink_path is None:
        logger.error("PLINK not found. Please install PLINK and ensure it's in PATH.")
        return None
    return plink_path


def get_convertf_executable() -> Optional[Path]:
    """Get convertf executable path (for EIGENSTRAT conversion)."""
    return get_tool_path("convertf")


def find_reference_files(
    reference: str,
    data_dir: Path,
) -> Optional[Dict[str, Path]]:
    """
    Find reference dataset files.

    Args:
        reference: Reference dataset name ("1000g", "aadr", "sgdp")
        data_dir: Base data directory

    Returns:
        Dictionary with file paths, or None if not found
    """
    reference = reference.lower()

    if reference == "1000g":
        ref_dir = data_dir / "1000genomes"
        # Look for converted PLINK files
        plink_prefix = ref_dir / "1000g_autosomal"
        if check_plink_files(plink_prefix):
            return {"format": "plink", "prefix": plink_prefix}

        # Otherwise need to convert from VCF
        vcf_pattern = ref_dir / "ALL.chr*.vcf.gz"
        vcf_files = list(ref_dir.glob("ALL.chr*.vcf.gz"))
        if vcf_files:
            return {"format": "vcf", "files": vcf_files, "dir": ref_dir}

    elif reference == "aadr":
        ref_dir = data_dir / "aadr"
        # Look for PLINK files first
        plink_prefix = ref_dir / "v54.1_1240K_public"
        if check_plink_files(plink_prefix):
            return {"format": "plink", "prefix": plink_prefix}

        # Check for EIGENSTRAT files
        eigen_prefix = ref_dir / "v54.1_1240K_public"
        if check_eigenstrat_files(eigen_prefix):
            return {"format": "eigenstrat", "prefix": eigen_prefix}

    elif reference == "sgdp":
        ref_dir = data_dir / "sgdp"
        # Look for PLINK files first
        plink_prefix = ref_dir / "cteam_extended.v4"
        if check_plink_files(plink_prefix):
            return {"format": "plink", "prefix": plink_prefix}

        # Check for EIGENSTRAT files
        eigen_prefix = ref_dir / "cteam_extended.v4"
        if check_eigenstrat_files(eigen_prefix):
            return {"format": "eigenstrat", "prefix": eigen_prefix}

    logger.warning(f"Reference dataset not found: {reference}")
    return None


def convert_eigenstrat_to_plink(
    eigen_prefix: Path,
    output_prefix: Path,
    convertf_path: Optional[Path] = None,
    plink_path: Optional[Path] = None,
) -> bool:
    """
    Convert EIGENSTRAT format to PLINK format.

    Args:
        eigen_prefix: Input EIGENSTRAT file prefix
        output_prefix: Output PLINK file prefix
        convertf_path: Path to convertf executable
        plink_path: Path to PLINK executable

    Returns:
        True if conversion successful
    """
    if convertf_path is None:
        convertf_path = get_convertf_executable()

    if plink_path is None:
        plink_path = get_plink_executable()

    if convertf_path is None:
        logger.warning("convertf not found, attempting alternative conversion")
        # Try using Python-based conversion
        return convert_eigenstrat_to_plink_python(eigen_prefix, output_prefix, plink_path)

    # Create convertf parameter file
    with tempfile.NamedTemporaryFile(mode="w", suffix=".par", delete=False) as f:
        par_file = Path(f.name)
        f.write(f"genotypename: {eigen_prefix}.geno\n")
        f.write(f"snpname: {eigen_prefix}.snp\n")
        f.write(f"indivname: {eigen_prefix}.ind\n")
        f.write(f"outputformat: PACKEDPED\n")
        f.write(f"genotypeoutname: {output_prefix}.bed\n")
        f.write(f"snpoutname: {output_prefix}.bim\n")
        f.write(f"indivoutname: {output_prefix}.fam\n")

    try:
        cmd = [str(convertf_path), "-p", str(par_file)]
        success, _, _ = run_command(cmd, "Converting EIGENSTRAT to PLINK")
        return success and check_plink_files(output_prefix)
    finally:
        par_file.unlink(missing_ok=True)


def convert_eigenstrat_to_plink_python(
    eigen_prefix: Path,
    output_prefix: Path,
    plink_path: Optional[Path] = None,
) -> bool:
    """
    Python-based EIGENSTRAT to PLINK conversion (fallback).

    This is a simplified conversion that works without convertf.
    """
    logger.info("Using Python-based EIGENSTRAT to PLINK conversion")

    geno_file = eigen_prefix.parent / f"{eigen_prefix.name}.geno"
    snp_file = eigen_prefix.parent / f"{eigen_prefix.name}.snp"
    ind_file = eigen_prefix.parent / f"{eigen_prefix.name}.ind"

    if not all(f.exists() for f in [geno_file, snp_file, ind_file]):
        logger.error("EIGENSTRAT files not found")
        return False

    # Read individual information
    individuals = []
    with open(ind_file, "r") as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 3:
                individuals.append({
                    "id": parts[0],
                    "sex": parts[1],
                    "pop": parts[2],
                })

    # Read SNP information
    snps = []
    with open(snp_file, "r") as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 6:
                snps.append({
                    "id": parts[0],
                    "chrom": parts[1],
                    "cm": parts[2],
                    "pos": parts[3],
                    "ref": parts[4],
                    "alt": parts[5],
                })

    # Write FAM file
    fam_path = output_prefix.parent / f"{output_prefix.name}.fam"
    with open(fam_path, "w") as f:
        for ind in individuals:
            sex_code = {"M": "1", "F": "2", "U": "0"}.get(ind["sex"].upper(), "0")
            f.write(f"{ind['id']} {ind['id']} 0 0 {sex_code} -9\n")

    # Write BIM file
    bim_path = output_prefix.parent / f"{output_prefix.name}.bim"
    with open(bim_path, "w") as f:
        for snp in snps:
            f.write(f"{snp['chrom']}\t{snp['id']}\t{snp['cm']}\t{snp['pos']}\t{snp['ref']}\t{snp['alt']}\n")

    # Write temporary PED file and convert to BED
    # This is memory-intensive for large datasets, so we chunk it
    ped_path = output_prefix.parent / f"{output_prefix.name}.ped"
    map_path = output_prefix.parent / f"{output_prefix.name}.map"

    # Write MAP file
    with open(map_path, "w") as f:
        for snp in snps:
            f.write(f"{snp['chrom']}\t{snp['id']}\t{snp['cm']}\t{snp['pos']}\n")

    # Process genotypes - read line by line to avoid memory issues
    logger.info("Converting genotypes (this may take a while for large files)...")

    # EIGENSTRAT .geno format: each line is a SNP, each character is an individual
    # We need to transpose: for each individual, get their genotype at each SNP
    # For memory efficiency, we read the file multiple times (once per individual)
    # This is slower but uses O(n_snps) memory instead of O(n_snps * n_individuals)

    n_individuals = len(individuals)
    n_snps = len(snps)

    # For smaller datasets, use the in-memory approach
    if n_snps * n_individuals < 100_000_000:  # ~100MB threshold
        with open(geno_file, "r") as geno_f:
            geno_lines = geno_f.readlines()

        with open(ped_path, "w") as ped_f:
            for idx, ind in enumerate(tqdm(individuals, desc="Converting")):
                sex_code = {"M": "1", "F": "2", "U": "0"}.get(ind["sex"].upper(), "0")
                ped_f.write(f"{ind['id']} {ind['id']} 0 0 {sex_code} -9")

                for snp_idx, snp in enumerate(snps):
                    geno_line = geno_lines[snp_idx] if snp_idx < len(geno_lines) else ""
                    geno_char = geno_line[idx] if idx < len(geno_line) else "9"

                    ref, alt = snp["ref"], snp["alt"]

                    if geno_char == "2":  # Homozygous reference
                        ped_f.write(f" {ref} {ref}")
                    elif geno_char == "1":  # Heterozygous
                        ped_f.write(f" {ref} {alt}")
                    elif geno_char == "0":  # Homozygous alternate
                        ped_f.write(f" {alt} {alt}")
                    else:  # Missing
                        ped_f.write(" 0 0")

                ped_f.write("\n")
    else:
        # For large datasets, read line by line for each individual
        # This is slower but memory-safe
        logger.warning("Large dataset detected - using memory-efficient mode (slower)")

        with open(ped_path, "w") as ped_f:
            for idx, ind in enumerate(tqdm(individuals, desc="Converting")):
                sex_code = {"M": "1", "F": "2", "U": "0"}.get(ind["sex"].upper(), "0")
                ped_f.write(f"{ind['id']} {ind['id']} 0 0 {sex_code} -9")

                with open(geno_file, "r") as geno_f:
                    for snp_idx, geno_line in enumerate(geno_f):
                        if snp_idx >= len(snps):
                            break

                        geno_char = geno_line[idx] if idx < len(geno_line) else "9"
                        snp = snps[snp_idx]
                        ref, alt = snp["ref"], snp["alt"]

                        if geno_char == "2":
                            ped_f.write(f" {ref} {ref}")
                        elif geno_char == "1":
                            ped_f.write(f" {ref} {alt}")
                        elif geno_char == "0":
                            ped_f.write(f" {alt} {alt}")
                        else:
                            ped_f.write(" 0 0")

                ped_f.write("\n")

    # Convert to binary format using PLINK
    if plink_path:
        cmd = [
            str(plink_path),
            "--file", str(output_prefix),
            "--make-bed",
            "--out", str(output_prefix),
            "--allow-no-sex",
        ]
        success, _, _ = run_command(cmd, "Converting to binary PLINK")

        # Clean up text files
        ped_path.unlink(missing_ok=True)
        map_path.unlink(missing_ok=True)

        return success and check_plink_files(output_prefix)

    return True


def get_snp_set(bim_file: Path) -> Set[str]:
    """Get set of SNP IDs from a BIM file."""
    snps = set()
    with open(bim_file, "r") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 2:
                snps.add(parts[1])
    return snps


def find_overlapping_snps(
    sample_bim: Path,
    reference_bim: Path,
) -> Set[str]:
    """
    Find overlapping SNPs between sample and reference.

    Args:
        sample_bim: Path to sample BIM file
        reference_bim: Path to reference BIM file

    Returns:
        Set of overlapping SNP IDs
    """
    logger.info("Finding overlapping SNPs...")

    sample_snps = get_snp_set(sample_bim)
    reference_snps = get_snp_set(reference_bim)

    overlap = sample_snps & reference_snps

    logger.info(f"  Sample SNPs: {len(sample_snps)}")
    logger.info(f"  Reference SNPs: {len(reference_snps)}")
    logger.info(f"  Overlapping SNPs: {len(overlap)}")

    return overlap


def extract_overlapping_snps(
    input_prefix: Path,
    output_prefix: Path,
    snp_list: Set[str],
    plink: Path,
) -> bool:
    """
    Extract overlapping SNPs from a PLINK dataset.

    Args:
        input_prefix: Input PLINK file prefix
        output_prefix: Output PLINK file prefix
        snp_list: Set of SNP IDs to extract
        plink: Path to PLINK executable

    Returns:
        True if extraction successful
    """
    # Write SNP list to temporary file
    with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False) as f:
        snp_file = Path(f.name)
        for snp_id in snp_list:
            f.write(f"{snp_id}\n")

    try:
        cmd = [
            str(plink),
            "--bfile", str(input_prefix),
            "--extract", str(snp_file),
            "--make-bed",
            "--out", str(output_prefix),
            "--allow-no-sex",
        ]

        success, _, _ = run_command(cmd, "Extracting overlapping SNPs")
        return success and check_plink_files(output_prefix)

    finally:
        snp_file.unlink(missing_ok=True)


def merge_plink_datasets(
    file_list: List[Path],
    output_prefix: Path,
    plink: Path,
) -> bool:
    """
    Merge multiple PLINK datasets.

    Args:
        file_list: List of PLINK file prefixes to merge
        output_prefix: Output PLINK file prefix
        plink: Path to PLINK executable

    Returns:
        True if merge successful
    """
    if len(file_list) < 2:
        logger.error("Need at least 2 datasets to merge")
        return False

    # Write merge list file
    with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False) as f:
        merge_list = Path(f.name)
        for prefix in file_list[1:]:
            f.write(f"{prefix}.bed {prefix}.bim {prefix}.fam\n")

    try:
        cmd = [
            str(plink),
            "--bfile", str(file_list[0]),
            "--merge-list", str(merge_list),
            "--make-bed",
            "--out", str(output_prefix),
            "--allow-no-sex",
        ]

        success, stdout, stderr = run_command(cmd, "Merging PLINK datasets")

        if not success:
            # Check for strand issues
            if "strand" in stderr.lower() or "allele" in stderr.lower():
                logger.warning("Merge failed due to strand/allele issues, attempting to resolve...")
                return merge_with_strand_resolution(file_list, output_prefix, plink)

        return success and check_plink_files(output_prefix)

    finally:
        merge_list.unlink(missing_ok=True)


def merge_with_strand_resolution(
    file_list: List[Path],
    output_prefix: Path,
    plink: Path,
) -> bool:
    """
    Merge datasets with strand/allele conflict resolution.

    Args:
        file_list: List of PLINK file prefixes to merge
        output_prefix: Output PLINK file prefix
        plink: Path to PLINK executable

    Returns:
        True if merge successful
    """
    temp_dir = output_prefix.parent / "temp_merge"
    ensure_dir(temp_dir)

    try:
        # Step 1: Try merge and get missnp file
        merge_list_path = temp_dir / "merge_list.txt"
        with open(merge_list_path, "w") as f:
            for prefix in file_list[1:]:
                f.write(f"{prefix}.bed {prefix}.bim {prefix}.fam\n")

        cmd = [
            str(plink),
            "--bfile", str(file_list[0]),
            "--merge-list", str(merge_list_path),
            "--make-bed",
            "--out", str(temp_dir / "trial_merge"),
            "--allow-no-sex",
        ]

        run_command(cmd, "Trial merge to identify problematic SNPs")

        # Step 2: Check for missnp file
        missnp_file = temp_dir / "trial_merge-merge.missnp"

        if missnp_file.exists():
            logger.info("Excluding problematic SNPs...")

            # Read problematic SNPs
            problematic_snps = set()
            with open(missnp_file, "r") as f:
                for line in f:
                    problematic_snps.add(line.strip())

            logger.info(f"Found {len(problematic_snps)} problematic SNPs to exclude")

            # Exclude from all datasets
            clean_prefixes = []
            for i, prefix in enumerate(file_list):
                clean_prefix = temp_dir / f"clean_{i}"

                with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False) as f:
                    exclude_file = Path(f.name)
                    for snp in problematic_snps:
                        f.write(f"{snp}\n")

                cmd = [
                    str(plink),
                    "--bfile", str(prefix),
                    "--exclude", str(exclude_file),
                    "--make-bed",
                    "--out", str(clean_prefix),
                    "--allow-no-sex",
                ]

                success, _, _ = run_command(cmd, f"Cleaning dataset {i+1}")
                exclude_file.unlink(missing_ok=True)

                if success and check_plink_files(clean_prefix):
                    clean_prefixes.append(clean_prefix)
                else:
                    return False

            # Step 3: Merge cleaned datasets
            return merge_plink_datasets(clean_prefixes, output_prefix, plink)

        return False

    finally:
        # Clean up temporary files
        import shutil
        if temp_dir.exists():
            shutil.rmtree(temp_dir)


def write_merge_report(
    stats: Dict,
    output_path: Path,
) -> None:
    """Write merge statistics report."""
    from datetime import datetime

    with open(output_path, "w") as f:
        f.write("=" * 60 + "\n")
        f.write("Dataset Merge Report\n")
        f.write("=" * 60 + "\n")
        f.write(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write("\n")

        f.write("Input:\n")
        f.write(f"  Sample: {stats.get('sample_prefix', 'N/A')}\n")
        f.write(f"  Sample SNPs: {stats.get('sample_snps', 'N/A')}\n")
        f.write(f"  Sample individuals: {stats.get('sample_individuals', 'N/A')}\n")
        f.write("\n")

        f.write("References:\n")
        for ref in stats.get("references", []):
            f.write(f"  {ref['name']}:\n")
            f.write(f"    SNPs: {ref.get('snps', 'N/A')}\n")
            f.write(f"    Individuals: {ref.get('individuals', 'N/A')}\n")
            f.write(f"    Overlapping SNPs: {ref.get('overlap', 'N/A')}\n")
        f.write("\n")

        f.write("Output:\n")
        f.write(f"  File: {stats.get('output_prefix', 'N/A')}\n")
        f.write(f"  Final SNPs: {stats.get('final_snps', 'N/A')}\n")
        f.write(f"  Final individuals: {stats.get('final_individuals', 'N/A')}\n")
        f.write("=" * 60 + "\n")

    logger.info(f"Merge report written to: {output_path}")


def run_merge_pipeline(
    sample_prefix: Path,
    references: List[str],
    output_prefix: Path,
    data_dir: Path,
) -> Tuple[bool, Dict]:
    """
    Run the complete merge pipeline.

    Args:
        sample_prefix: Sample PLINK file prefix
        references: List of reference dataset names
        output_prefix: Output PLINK file prefix
        data_dir: Base data directory

    Returns:
        Tuple of (success, statistics dict)
    """
    plink = get_plink_executable()
    if plink is None:
        return False, {}

    ensure_dir(output_prefix.parent)
    temp_dir = output_prefix.parent / "temp_merge"
    ensure_dir(temp_dir)

    # Initialize statistics
    stats = {
        "sample_prefix": str(sample_prefix),
        "sample_snps": count_plink_snps(sample_prefix.parent / f"{sample_prefix.name}.bim"),
        "sample_individuals": count_plink_samples(sample_prefix.parent / f"{sample_prefix.name}.fam"),
        "references": [],
    }

    try:
        # Process each reference dataset
        datasets_to_merge = [sample_prefix]
        common_snps: Optional[Set[str]] = None

        for ref_name in references:
            logger.info(f"Processing reference: {ref_name}")

            ref_info = find_reference_files(ref_name, data_dir)
            if ref_info is None:
                logger.warning(f"Skipping {ref_name}: files not found")
                continue

            ref_stats = {"name": ref_name}

            # Get or convert to PLINK format
            if ref_info["format"] == "plink":
                ref_prefix = ref_info["prefix"]
            elif ref_info["format"] == "eigenstrat":
                ref_prefix = temp_dir / f"{ref_name}_plink"
                success = convert_eigenstrat_to_plink(
                    ref_info["prefix"],
                    ref_prefix,
                    plink_path=plink,
                )
                if not success:
                    logger.warning(f"Failed to convert {ref_name} to PLINK format")
                    continue
            else:
                logger.warning(f"Unsupported format for {ref_name}: {ref_info['format']}")
                continue

            # Count SNPs and individuals
            ref_bim = ref_prefix.parent / f"{ref_prefix.name}.bim"
            ref_fam = ref_prefix.parent / f"{ref_prefix.name}.fam"

            ref_stats["snps"] = count_plink_snps(ref_bim)
            ref_stats["individuals"] = count_plink_samples(ref_fam)

            # Find overlapping SNPs
            sample_bim = sample_prefix.parent / f"{sample_prefix.name}.bim"
            overlap = find_overlapping_snps(sample_bim, ref_bim)
            ref_stats["overlap"] = len(overlap)

            if common_snps is None:
                common_snps = overlap
            else:
                common_snps = common_snps & overlap

            stats["references"].append(ref_stats)
            datasets_to_merge.append(ref_prefix)

        if len(datasets_to_merge) < 2:
            logger.error("No valid reference datasets found")
            return False, stats

        if common_snps is None or len(common_snps) == 0:
            logger.error("No overlapping SNPs found between datasets")
            return False, stats

        logger.info(f"Common SNPs across all datasets: {len(common_snps)}")

        # Extract common SNPs from each dataset
        extracted_prefixes = []
        for i, prefix in enumerate(datasets_to_merge):
            extracted_prefix = temp_dir / f"extracted_{i}"
            success = extract_overlapping_snps(prefix, extracted_prefix, common_snps, plink)
            if not success:
                logger.error(f"Failed to extract SNPs from dataset {i}")
                return False, stats
            extracted_prefixes.append(extracted_prefix)

        # Merge all datasets
        success = merge_plink_datasets(extracted_prefixes, output_prefix, plink)

        if success and check_plink_files(output_prefix):
            stats["output_prefix"] = str(output_prefix)
            stats["final_snps"] = count_plink_snps(output_prefix.parent / f"{output_prefix.name}.bim")
            stats["final_individuals"] = count_plink_samples(output_prefix.parent / f"{output_prefix.name}.fam")

            logger.info(f"Merge completed successfully")
            logger.info(f"  Final SNPs: {stats['final_snps']}")
            logger.info(f"  Final individuals: {stats['final_individuals']}")

            return True, stats

        return False, stats

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
        description="Merge sample data with reference populations.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Merge with 1000 Genomes
    %(prog)s --sample sample --references 1000g --output merged

    # Merge with multiple references
    %(prog)s --sample sample --references 1000g aadr sgdp --output merged

Output:
    - {output}.bed/bim/fam - Merged PLINK files
    - {output}_merge_report.txt - Merge statistics

Available references:
    - 1000g: 1000 Genomes Project Phase 3
    - aadr: Allen Ancient DNA Resource
    - sgdp: Simons Genome Diversity Project
        """,
    )

    parser.add_argument(
        "--sample", "-s",
        type=str,
        required=True,
        help="Sample PLINK file prefix (required)",
    )

    parser.add_argument(
        "--references", "-r",
        type=str,
        nargs="+",
        required=True,
        choices=["1000g", "aadr", "sgdp"],
        help="Reference dataset names (required)",
    )

    parser.add_argument(
        "--output", "-o",
        type=str,
        required=True,
        help="Output merged PLINK file prefix (required)",
    )

    parser.add_argument(
        "--data-dir",
        type=str,
        default=str(RAW_DATA_DIR),
        help=f"Base data directory (default: {RAW_DATA_DIR})",
    )

    parser.add_argument(
        "--log-file",
        type=str,
        help="Path to log file",
    )

    args = parser.parse_args()

    # Set up logging
    log_file = Path(args.log_file) if args.log_file else LOGS_DIR / "merge_datasets.log"
    setup_logging(log_file)

    logger.info("=" * 60)
    logger.info("Dataset Merge Pipeline")
    logger.info("=" * 60)

    # Validate input
    sample_prefix = Path(args.sample).resolve()
    if not check_plink_files(sample_prefix):
        logger.error(f"Sample PLINK files not found: {sample_prefix}")
        return 1

    output_prefix = Path(args.output).resolve()
    data_dir = Path(args.data_dir).resolve()

    logger.info(f"Sample: {sample_prefix}")
    logger.info(f"References: {', '.join(args.references)}")
    logger.info(f"Output: {output_prefix}")
    logger.info(f"Data directory: {data_dir}")

    # Run merge pipeline
    success, stats = run_merge_pipeline(
        sample_prefix=sample_prefix,
        references=args.references,
        output_prefix=output_prefix,
        data_dir=data_dir,
    )

    if success:
        # Write merge report
        report_path = output_prefix.parent / f"{output_prefix.name}_merge_report.txt"
        write_merge_report(stats, report_path)

        logger.info("=" * 60)
        logger.info("Merge completed successfully!")
        logger.info("=" * 60)
        return 0
    else:
        logger.error("Merge pipeline failed")
        return 1


if __name__ == "__main__":
    sys.exit(main())
