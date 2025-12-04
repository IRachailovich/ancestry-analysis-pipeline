#!/usr/bin/env python3
"""
Convert 23andMe raw data to PLINK format.

This script converts 23andMe genotype data from their text file format
to PLINK binary format (.bed/.bim/.fam).

Features:
    - Parses 23andMe text file format (rsID, chromosome, position, genotype)
    - Handles different 23andMe file versions (v3, v4, v5) - detects automatically
    - Converts to PLINK .ped/.map format, then to binary .bed/.bim/.fam
    - Handles missing genotypes (--,--, DD), non-autosomal SNPs (X, Y, MT)
    - Filters to only autosomal chromosomes (1-22)
    - Progress bars with tqdm
    - Comprehensive logging

Usage:
    python convert_23andme.py --input data/raw/23andme/genome.txt \\
                              --output data/processed/sample \\
                              --sample-id MYSAMPLE \\
                              --sex 2

Example:
    # Basic conversion
    python convert_23andme.py --input genome.txt --output sample

    # With custom sample ID and sex
    python convert_23andme.py --input genome.txt --output sample \\
                              --sample-id "SAMPLE001" --sex 1
"""

import argparse
import logging
import os
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from tqdm import tqdm

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from utils.config import LOGS_DIR, get_tool_path
from utils.file_utils import ensure_dir, check_plink_files
from utils.genetics_utils import (
    detect_23andme_version,
    is_autosomal,
    normalize_chromosome,
    parse_23andme_genotype,
    is_valid_rsid,
)


logger = logging.getLogger(__name__)


def setup_logging(log_file: Optional[Path] = None) -> logging.Logger:
    """
    Set up logging configuration.

    Args:
        log_file: Optional path to log file

    Returns:
        Configured logger
    """
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


def parse_23andme_file(
    input_path: Path,
    filter_autosomal: bool = True,
    show_progress: bool = True,
) -> Tuple[List[Dict], Dict]:
    """
    Parse a 23andMe raw data file.

    Args:
        input_path: Path to 23andMe raw data file
        filter_autosomal: If True, only keep autosomal chromosomes (1-22)
        show_progress: Whether to show progress bar

    Returns:
        Tuple of (list of SNP records, file statistics dict)

    Raises:
        FileNotFoundError: If input file doesn't exist
        ValueError: If file format is invalid
    """
    if not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_path}")

    # Read header lines to detect version
    header_lines = []
    data_lines = []

    logger.info(f"Reading 23andMe file: {input_path}")

    # First pass: count lines for progress bar
    total_lines = sum(1 for _ in open(input_path, "r", encoding="utf-8", errors="replace"))

    with open(input_path, "r", encoding="utf-8", errors="replace") as f:
        iterator = tqdm(f, total=total_lines, desc="Parsing", disable=not show_progress)

        for line in iterator:
            line = line.strip()

            # Skip empty lines
            if not line:
                continue

            # Collect header lines (comments)
            if line.startswith("#"):
                header_lines.append(line)
                continue

            data_lines.append(line)

    # Detect file version
    version = detect_23andme_version(header_lines)
    logger.info(f"Detected 23andMe file version: {version}")

    # Parse data lines
    snps = []
    stats = {
        "total_lines": len(data_lines),
        "valid_snps": 0,
        "missing_genotypes": 0,
        "non_autosomal": 0,
        "invalid_format": 0,
        "version": version,
    }

    for line in tqdm(data_lines, desc="Processing SNPs", disable=not show_progress):
        parts = line.split("\t")

        if len(parts) < 4:
            stats["invalid_format"] += 1
            continue

        rsid = parts[0].strip()
        chromosome = parts[1].strip()
        position = parts[2].strip()
        genotype = parts[3].strip()

        # Validate position
        try:
            pos_int = int(position)
            if pos_int <= 0:
                stats["invalid_format"] += 1
                continue
        except ValueError:
            stats["invalid_format"] += 1
            continue

        # Check chromosome
        if filter_autosomal:
            if not is_autosomal(chromosome):
                stats["non_autosomal"] += 1
                continue

        # Parse genotype
        allele1, allele2 = parse_23andme_genotype(genotype)

        # Check for missing genotypes
        if allele1 == "0" or allele2 == "0":
            stats["missing_genotypes"] += 1
            continue

        # Create SNP record
        snp_record = {
            "rsid": rsid,
            "chromosome": normalize_chromosome(chromosome),
            "position": pos_int,
            "allele1": allele1,
            "allele2": allele2,
        }

        snps.append(snp_record)
        stats["valid_snps"] += 1

    logger.info(f"Parsing complete:")
    logger.info(f"  Total lines: {stats['total_lines']}")
    logger.info(f"  Valid SNPs: {stats['valid_snps']}")
    logger.info(f"  Missing genotypes: {stats['missing_genotypes']}")
    logger.info(f"  Non-autosomal: {stats['non_autosomal']}")
    logger.info(f"  Invalid format: {stats['invalid_format']}")

    return snps, stats


def write_plink_ped_map(
    snps: List[Dict],
    output_prefix: Path,
    sample_id: str = "SAMPLE001",
    family_id: Optional[str] = None,
    sex: int = 0,
    phenotype: str = "-9",
    show_progress: bool = True,
) -> Tuple[Path, Path]:
    """
    Write PLINK .ped and .map files.

    Args:
        snps: List of SNP records from parse_23andme_file
        output_prefix: Output file prefix (without extension)
        sample_id: Sample identifier
        family_id: Family identifier (defaults to sample_id)
        sex: Sex code (1=male, 2=female, 0=unknown)
        phenotype: Phenotype value
        show_progress: Whether to show progress bar

    Returns:
        Tuple of (ped_path, map_path)
    """
    if family_id is None:
        family_id = sample_id

    output_prefix = Path(output_prefix)
    output_prefix.parent.mkdir(parents=True, exist_ok=True)

    ped_path = output_prefix.parent / f"{output_prefix.name}.ped"
    map_path = output_prefix.parent / f"{output_prefix.name}.map"

    logger.info(f"Writing PLINK .map file: {map_path}")

    # Write .map file
    with open(map_path, "w") as f:
        for snp in tqdm(snps, desc="Writing .map", disable=not show_progress):
            # Format: CHR SNP_ID CM BP
            f.write(f"{snp['chromosome']}\t{snp['rsid']}\t0\t{snp['position']}\n")

    logger.info(f"Writing PLINK .ped file: {ped_path}")

    # Write .ped file
    with open(ped_path, "w") as f:
        # First 6 columns: FID IID PAT MAT SEX PHENO
        f.write(f"{family_id} {sample_id} 0 0 {sex} {phenotype}")

        # Remaining columns: genotypes (space-separated allele pairs)
        for snp in tqdm(snps, desc="Writing .ped", disable=not show_progress):
            f.write(f" {snp['allele1']} {snp['allele2']}")
        f.write("\n")

    logger.info(f"Written {len(snps)} SNPs to PLINK files")

    return ped_path, map_path


def convert_to_binary_plink(
    ped_path: Path,
    output_prefix: Path,
    plink_path: Optional[Path] = None,
) -> bool:
    """
    Convert PLINK text format (.ped/.map) to binary format (.bed/.bim/.fam).

    Args:
        ped_path: Path to .ped file
        output_prefix: Output prefix for binary files
        plink_path: Path to PLINK executable (auto-detected if None)

    Returns:
        True if conversion successful, False otherwise
    """
    # Find PLINK executable
    if plink_path is None:
        plink_path = get_tool_path("plink")
        if plink_path is None:
            logger.warning("PLINK not found. Binary conversion skipped.")
            logger.info("Install PLINK and ensure it's in PATH to enable binary conversion.")
            return False

    map_path = ped_path.with_suffix(".map")
    output_prefix = Path(output_prefix)

    logger.info(f"Converting to binary PLINK format using: {plink_path}")

    # Build PLINK command
    cmd = [
        str(plink_path),
        "--ped", str(ped_path),
        "--map", str(map_path),
        "--make-bed",
        "--out", str(output_prefix),
        "--allow-no-sex",
    ]

    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=False,
        )

        if result.returncode != 0:
            logger.error(f"PLINK conversion failed: {result.stderr}")
            return False

        # Check output files exist
        if check_plink_files(output_prefix):
            logger.info("Binary PLINK files created successfully")

            # Clean up temporary files
            ped_path.unlink(missing_ok=True)
            map_path.unlink(missing_ok=True)

            # Also clean up PLINK log/nosex files
            log_file = output_prefix.parent / f"{output_prefix.name}.log"
            nosex_file = output_prefix.parent / f"{output_prefix.name}.nosex"
            log_file.unlink(missing_ok=True)
            nosex_file.unlink(missing_ok=True)

            return True
        else:
            logger.error("PLINK conversion completed but output files not found")
            return False

    except FileNotFoundError:
        logger.error(f"PLINK executable not found at: {plink_path}")
        return False
    except Exception as e:
        logger.error(f"Error running PLINK: {e}")
        return False


def convert_23andme_to_plink(
    input_path: Path,
    output_prefix: Path,
    sample_id: str = "SAMPLE001",
    sex: int = 0,
    convert_binary: bool = True,
    show_progress: bool = True,
) -> Tuple[bool, Dict]:
    """
    Convert 23andMe raw data to PLINK format.

    This is the main conversion function that orchestrates the full pipeline.

    Args:
        input_path: Path to 23andMe raw data file
        output_prefix: Output prefix for PLINK files
        sample_id: Sample identifier
        sex: Sex code (1=male, 2=female, 0=unknown)
        convert_binary: If True, convert to binary .bed/.bim/.fam format
        show_progress: Whether to show progress bars

    Returns:
        Tuple of (success, statistics_dict)
    """
    # Parse 23andMe file
    snps, stats = parse_23andme_file(
        input_path,
        filter_autosomal=True,
        show_progress=show_progress,
    )

    if not snps:
        logger.error("No valid SNPs found in input file")
        return False, stats

    # Write PLINK text format
    ped_path, map_path = write_plink_ped_map(
        snps=snps,
        output_prefix=output_prefix,
        sample_id=sample_id,
        sex=sex,
        show_progress=show_progress,
    )

    stats["ped_path"] = str(ped_path)
    stats["map_path"] = str(map_path)

    # Convert to binary format if requested
    if convert_binary:
        success = convert_to_binary_plink(ped_path, output_prefix)
        stats["binary_conversion"] = success

        if success:
            stats["bed_path"] = str(output_prefix.parent / f"{output_prefix.name}.bed")
            stats["bim_path"] = str(output_prefix.parent / f"{output_prefix.name}.bim")
            stats["fam_path"] = str(output_prefix.parent / f"{output_prefix.name}.fam")
    else:
        stats["binary_conversion"] = False

    return True, stats


def main() -> int:
    """
    Main entry point for the 23andMe conversion script.

    Returns:
        Exit code (0 for success, non-zero for failure)
    """
    parser = argparse.ArgumentParser(
        description="Convert 23andMe raw data to PLINK format.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Basic conversion
    %(prog)s --input genome.txt --output data/processed/sample

    # With custom sample ID and sex
    %(prog)s --input genome.txt --output sample --sample-id MYSAMPLE --sex 2

    # Keep text format (no binary conversion)
    %(prog)s --input genome.txt --output sample --no-binary

Output files:
    - {output}.bed - Binary genotype data
    - {output}.bim - SNP information
    - {output}.fam - Sample information
        """,
    )

    parser.add_argument(
        "--input", "-i",
        type=str,
        required=True,
        help="Path to 23andMe raw data file",
    )

    parser.add_argument(
        "--output", "-o",
        type=str,
        required=True,
        help="Output prefix for PLINK files",
    )

    parser.add_argument(
        "--sample-id",
        type=str,
        default="SAMPLE001",
        help="Sample ID (default: SAMPLE001)",
    )

    parser.add_argument(
        "--sex",
        type=int,
        choices=[0, 1, 2],
        default=0,
        help="Sex code: 1=male, 2=female, 0=unknown (default: 0)",
    )

    parser.add_argument(
        "--no-binary",
        action="store_true",
        help="Skip binary conversion (keep .ped/.map format)",
    )

    parser.add_argument(
        "--log-file",
        type=str,
        help="Path to log file (default: logs/convert_23andme.log)",
    )

    parser.add_argument(
        "--quiet", "-q",
        action="store_true",
        help="Suppress progress bars",
    )

    args = parser.parse_args()

    # Set up logging
    log_file = Path(args.log_file) if args.log_file else LOGS_DIR / "convert_23andme.log"
    setup_logging(log_file)

    logger.info("=" * 60)
    logger.info("23andMe to PLINK Converter")
    logger.info("=" * 60)

    # Validate input
    input_path = Path(args.input).resolve()
    if not input_path.exists():
        logger.error(f"Input file not found: {input_path}")
        return 1

    output_prefix = Path(args.output).resolve()

    logger.info(f"Input file: {input_path}")
    logger.info(f"Output prefix: {output_prefix}")
    logger.info(f"Sample ID: {args.sample_id}")
    logger.info(f"Sex: {args.sex}")
    logger.info(f"Binary conversion: {not args.no_binary}")

    # Run conversion
    try:
        success, stats = convert_23andme_to_plink(
            input_path=input_path,
            output_prefix=output_prefix,
            sample_id=args.sample_id,
            sex=args.sex,
            convert_binary=not args.no_binary,
            show_progress=not args.quiet,
        )

        if success:
            logger.info("=" * 60)
            logger.info("Conversion completed successfully!")
            logger.info(f"  Total SNPs processed: {stats['valid_snps']}")
            if stats.get("binary_conversion"):
                logger.info(f"  Output files: {stats.get('bed_path', 'N/A')}")
            else:
                logger.info(f"  Output files: {stats.get('ped_path', 'N/A')}")
            logger.info("=" * 60)
            return 0
        else:
            logger.error("Conversion failed")
            return 1

    except Exception as e:
        logger.exception(f"Unexpected error: {e}")
        return 1


if __name__ == "__main__":
    sys.exit(main())