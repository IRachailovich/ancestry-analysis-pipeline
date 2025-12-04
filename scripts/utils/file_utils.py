#!/usr/bin/env python3
"""
File handling utilities for genetic data formats.

This module provides utilities for:
- Checking if PLINK files exist (.bed, .bim, .fam)
- Checking if EIGENSTRAT files exist (.geno, .snp, .ind)
- Creating directory structures
- Reading/writing genetic data files
- Progress bars for file operations
- File size formatting

Example:
    >>> from scripts.utils.file_utils import check_plink_files, ensure_dir
    >>> if check_plink_files("data/processed/sample"):
    ...     print("PLINK files found")
    >>> ensure_dir("results/pca")
"""

import logging
import os
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

import pandas as pd
from tqdm import tqdm


logger = logging.getLogger(__name__)


# ==============================================================================
# Directory Utilities
# ==============================================================================

def ensure_dir(path: Union[str, Path]) -> Path:
    """
    Ensure a directory exists, creating it if necessary.

    Args:
        path: Path to the directory

    Returns:
        Path object for the directory

    Example:
        >>> output_dir = ensure_dir("results/pca")
        >>> output_dir.exists()
        True
    """
    dir_path = Path(path)
    dir_path.mkdir(parents=True, exist_ok=True)
    return dir_path


def create_project_structure(base_dir: Union[str, Path]) -> Dict[str, Path]:
    """
    Create the full project directory structure.

    Args:
        base_dir: Base project directory

    Returns:
        Dictionary mapping directory names to their paths
    """
    base = Path(base_dir)

    directories = {
        "data": base / "data",
        "raw": base / "data" / "raw",
        "raw_23andme": base / "data" / "raw" / "23andme",
        "raw_1000g": base / "data" / "raw" / "1000genomes",
        "raw_aadr": base / "data" / "raw" / "aadr",
        "raw_sgdp": base / "data" / "raw" / "sgdp",
        "processed": base / "data" / "processed",
        "reference_panels": base / "data" / "reference_panels",
        "results": base / "results",
        "pca": base / "results" / "pca",
        "admixture": base / "results" / "admixture",
        "qpadm": base / "results" / "qpadm",
        "local_ancestry": base / "results" / "local_ancestry",
        "admixture_graphs": base / "results" / "admixture_graphs",
        "figures": base / "results" / "figures",
        "logs": base / "logs",
    }

    for name, path in directories.items():
        path.mkdir(parents=True, exist_ok=True)

    return directories


# ==============================================================================
# PLINK File Utilities
# ==============================================================================

def check_plink_files(prefix: Union[str, Path]) -> bool:
    """
    Check if all PLINK binary files exist.

    Args:
        prefix: Path prefix for PLINK files (without extension)

    Returns:
        True if all files (.bed, .bim, .fam) exist, False otherwise

    Example:
        >>> check_plink_files("data/processed/sample")
        True
    """
    prefix = Path(prefix)
    extensions = [".bed", ".bim", ".fam"]
    return all((prefix.parent / f"{prefix.name}{ext}").exists() for ext in extensions)


def get_plink_file_paths(prefix: Union[str, Path]) -> Dict[str, Path]:
    """
    Get paths to all PLINK files for a given prefix.

    Args:
        prefix: Path prefix for PLINK files

    Returns:
        Dictionary with 'bed', 'bim', 'fam' keys and Path values
    """
    prefix = Path(prefix)
    return {
        "bed": prefix.parent / f"{prefix.name}.bed",
        "bim": prefix.parent / f"{prefix.name}.bim",
        "fam": prefix.parent / f"{prefix.name}.fam",
    }


def count_plink_snps(bim_file: Union[str, Path]) -> int:
    """
    Count the number of SNPs in a PLINK .bim file.

    Args:
        bim_file: Path to .bim file

    Returns:
        Number of SNPs
    """
    bim_path = Path(bim_file)
    if not bim_path.exists():
        raise FileNotFoundError(f"BIM file not found: {bim_path}")

    with open(bim_path, "r") as f:
        return sum(1 for _ in f)


def count_plink_samples(fam_file: Union[str, Path]) -> int:
    """
    Count the number of samples in a PLINK .fam file.

    Args:
        fam_file: Path to .fam file

    Returns:
        Number of samples
    """
    fam_path = Path(fam_file)
    if not fam_path.exists():
        raise FileNotFoundError(f"FAM file not found: {fam_path}")

    with open(fam_path, "r") as f:
        return sum(1 for _ in f)


def read_bim_file(bim_file: Union[str, Path]) -> pd.DataFrame:
    """
    Read a PLINK .bim file into a DataFrame.

    Args:
        bim_file: Path to .bim file

    Returns:
        DataFrame with columns: CHR, SNP, CM, BP, A1, A2
    """
    columns = ["CHR", "SNP", "CM", "BP", "A1", "A2"]
    return pd.read_csv(
        bim_file,
        sep="\t",
        header=None,
        names=columns,
        dtype={"CHR": str, "SNP": str, "CM": float, "BP": int, "A1": str, "A2": str},
    )


def read_fam_file(fam_file: Union[str, Path]) -> pd.DataFrame:
    """
    Read a PLINK .fam file into a DataFrame.

    Args:
        fam_file: Path to .fam file

    Returns:
        DataFrame with columns: FID, IID, PAT, MAT, SEX, PHENO
    """
    columns = ["FID", "IID", "PAT", "MAT", "SEX", "PHENO"]
    return pd.read_csv(
        fam_file,
        sep=r"\s+",
        header=None,
        names=columns,
        dtype={"FID": str, "IID": str, "PAT": str, "MAT": str, "SEX": int, "PHENO": str},
    )


def write_bim_file(df: pd.DataFrame, bim_file: Union[str, Path]) -> None:
    """
    Write a DataFrame to a PLINK .bim file.

    Args:
        df: DataFrame with BIM columns
        bim_file: Output path
    """
    df.to_csv(bim_file, sep="\t", header=False, index=False)


def write_fam_file(df: pd.DataFrame, fam_file: Union[str, Path]) -> None:
    """
    Write a DataFrame to a PLINK .fam file.

    Args:
        df: DataFrame with FAM columns
        fam_file: Output path
    """
    df.to_csv(fam_file, sep=" ", header=False, index=False)


def write_map_file(
    snp_ids: List[str],
    chromosomes: List[str],
    positions: List[int],
    output_path: Union[str, Path],
    cm_positions: Optional[List[float]] = None,
) -> None:
    """
    Write a PLINK .map file.

    Args:
        snp_ids: List of SNP IDs
        chromosomes: List of chromosome values
        positions: List of base pair positions
        output_path: Output file path
        cm_positions: Optional centiMorgan positions (defaults to 0)
    """
    if cm_positions is None:
        cm_positions = [0.0] * len(snp_ids)

    with open(output_path, "w") as f:
        for chrom, snp_id, cm, bp in zip(chromosomes, snp_ids, cm_positions, positions):
            f.write(f"{chrom}\t{snp_id}\t{cm}\t{bp}\n")


def write_ped_file(
    sample_id: str,
    family_id: str,
    sex: int,
    genotypes: List[Tuple[str, str]],
    output_path: Union[str, Path],
    paternal_id: str = "0",
    maternal_id: str = "0",
    phenotype: str = "-9",
) -> None:
    """
    Write a PLINK .ped file for a single sample.

    Args:
        sample_id: Sample ID
        family_id: Family ID
        sex: Sex code (1=male, 2=female, 0=unknown)
        genotypes: List of (allele1, allele2) tuples
        output_path: Output file path
        paternal_id: Paternal ID (default: "0")
        maternal_id: Maternal ID (default: "0")
        phenotype: Phenotype value (default: "-9")
    """
    with open(output_path, "w") as f:
        # First 6 columns: FID IID PAT MAT SEX PHENO
        f.write(f"{family_id} {sample_id} {paternal_id} {maternal_id} {sex} {phenotype}")

        # Remaining columns: genotypes (space-separated allele pairs)
        for a1, a2 in genotypes:
            f.write(f" {a1} {a2}")
        f.write("\n")


# ==============================================================================
# EIGENSTRAT File Utilities
# ==============================================================================

def check_eigenstrat_files(prefix: Union[str, Path]) -> bool:
    """
    Check if all EIGENSTRAT files exist.

    Args:
        prefix: Path prefix for EIGENSTRAT files (without extension)

    Returns:
        True if all files (.geno, .snp, .ind) exist, False otherwise

    Example:
        >>> check_eigenstrat_files("data/raw/aadr/v54.1_1240K_public")
        True
    """
    prefix = Path(prefix)
    extensions = [".geno", ".snp", ".ind"]
    return all((prefix.parent / f"{prefix.name}{ext}").exists() for ext in extensions)


def get_eigenstrat_file_paths(prefix: Union[str, Path]) -> Dict[str, Path]:
    """
    Get paths to all EIGENSTRAT files for a given prefix.

    Args:
        prefix: Path prefix for EIGENSTRAT files

    Returns:
        Dictionary with 'geno', 'snp', 'ind' keys and Path values
    """
    prefix = Path(prefix)
    return {
        "geno": prefix.parent / f"{prefix.name}.geno",
        "snp": prefix.parent / f"{prefix.name}.snp",
        "ind": prefix.parent / f"{prefix.name}.ind",
    }


def read_ind_file(ind_file: Union[str, Path]) -> pd.DataFrame:
    """
    Read an EIGENSTRAT .ind file into a DataFrame.

    Args:
        ind_file: Path to .ind file

    Returns:
        DataFrame with columns: ID, SEX, POPULATION
    """
    return pd.read_csv(
        ind_file,
        sep=r"\s+",
        header=None,
        names=["ID", "SEX", "POPULATION"],
        dtype=str,
    )


def read_snp_file(snp_file: Union[str, Path]) -> pd.DataFrame:
    """
    Read an EIGENSTRAT .snp file into a DataFrame.

    Args:
        snp_file: Path to .snp file

    Returns:
        DataFrame with columns: SNP, CHR, CM, BP, REF, ALT
    """
    return pd.read_csv(
        snp_file,
        sep=r"\s+",
        header=None,
        names=["SNP", "CHR", "CM", "BP", "REF", "ALT"],
        dtype={"SNP": str, "CHR": str, "CM": float, "BP": int, "REF": str, "ALT": str},
    )


def write_ind_file(
    sample_ids: List[str],
    sexes: List[str],
    populations: List[str],
    output_path: Union[str, Path],
) -> None:
    """
    Write an EIGENSTRAT .ind file.

    Args:
        sample_ids: List of sample IDs
        sexes: List of sex codes (M, F, U)
        populations: List of population labels
        output_path: Output file path
    """
    with open(output_path, "w") as f:
        for sid, sex, pop in zip(sample_ids, sexes, populations):
            f.write(f"{sid}\t{sex}\t{pop}\n")


def write_snp_file(
    snp_ids: List[str],
    chromosomes: List[str],
    positions: List[int],
    ref_alleles: List[str],
    alt_alleles: List[str],
    output_path: Union[str, Path],
    cm_positions: Optional[List[float]] = None,
) -> None:
    """
    Write an EIGENSTRAT .snp file.

    Args:
        snp_ids: List of SNP IDs
        chromosomes: List of chromosome values
        positions: List of base pair positions
        ref_alleles: List of reference alleles
        alt_alleles: List of alternate alleles
        output_path: Output file path
        cm_positions: Optional centiMorgan positions (defaults to 0)
    """
    if cm_positions is None:
        cm_positions = [0.0] * len(snp_ids)

    with open(output_path, "w") as f:
        for snp_id, chrom, cm, bp, ref, alt in zip(
            snp_ids, chromosomes, cm_positions, positions, ref_alleles, alt_alleles
        ):
            f.write(f"{snp_id}\t{chrom}\t{cm:.6f}\t{bp}\t{ref}\t{alt}\n")


# ==============================================================================
# File Size and Progress Utilities
# ==============================================================================

def format_size(size_bytes: int) -> str:
    """
    Format file size in human-readable format.

    Args:
        size_bytes: Size in bytes

    Returns:
        Human-readable size string (e.g., "1.23 GB")
    """
    for unit in ["B", "KB", "MB", "GB", "TB"]:
        if size_bytes < 1024.0:
            return f"{size_bytes:.2f} {unit}"
        size_bytes /= 1024.0
    return f"{size_bytes:.2f} PB"


def get_file_size(path: Union[str, Path]) -> int:
    """
    Get the size of a file in bytes.

    Args:
        path: Path to the file

    Returns:
        File size in bytes
    """
    return Path(path).stat().st_size


def count_lines(file_path: Union[str, Path], show_progress: bool = False) -> int:
    """
    Count the number of lines in a file.

    Args:
        file_path: Path to the file
        show_progress: Whether to show a progress bar

    Returns:
        Number of lines
    """
    file_path = Path(file_path)

    if show_progress:
        file_size = file_path.stat().st_size
        with open(file_path, "r") as f:
            with tqdm(total=file_size, unit="B", unit_scale=True, desc="Counting lines") as pbar:
                count = 0
                for line in f:
                    count += 1
                    pbar.update(len(line.encode("utf-8")))
        return count
    else:
        with open(file_path, "r") as f:
            return sum(1 for _ in f)


def read_file_with_progress(
    file_path: Union[str, Path],
    skip_comments: bool = True,
    comment_char: str = "#",
) -> List[str]:
    """
    Read a file line by line with a progress bar.

    Args:
        file_path: Path to the file
        skip_comments: Whether to skip comment lines
        comment_char: Character indicating comment lines

    Returns:
        List of lines (excluding comments if specified)
    """
    file_path = Path(file_path)
    file_size = file_path.stat().st_size

    lines = []
    with open(file_path, "r") as f:
        with tqdm(total=file_size, unit="B", unit_scale=True, desc=file_path.name) as pbar:
            for line in f:
                pbar.update(len(line.encode("utf-8")))
                if skip_comments and line.startswith(comment_char):
                    continue
                lines.append(line.rstrip("\n"))

    return lines


# ==============================================================================
# Population File Utilities
# ==============================================================================

def read_population_file(pop_file: Union[str, Path]) -> Dict[str, str]:
    """
    Read a population mapping file.

    Expected format: TAB or space separated, columns: SAMPLE_ID POPULATION
    Lines starting with # are treated as comments.

    Args:
        pop_file: Path to population file

    Returns:
        Dictionary mapping sample IDs to population labels
    """
    pop_map = {}
    with open(pop_file, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) >= 2:
                pop_map[parts[0]] = parts[1]
    return pop_map


def write_population_file(
    pop_map: Dict[str, str],
    output_path: Union[str, Path],
) -> None:
    """
    Write a population mapping file.

    Args:
        pop_map: Dictionary mapping sample IDs to population labels
        output_path: Output file path
    """
    with open(output_path, "w") as f:
        f.write("# Sample ID\tPopulation\n")
        for sample_id, population in sorted(pop_map.items()):
            f.write(f"{sample_id}\t{population}\n")


# ==============================================================================
# Extract/Filter Utilities
# ==============================================================================

def write_snp_list(snp_ids: List[str], output_path: Union[str, Path]) -> None:
    """
    Write a list of SNP IDs to a file.

    Args:
        snp_ids: List of SNP IDs
        output_path: Output file path
    """
    with open(output_path, "w") as f:
        for snp_id in snp_ids:
            f.write(f"{snp_id}\n")


def write_sample_list(
    sample_ids: List[str],
    output_path: Union[str, Path],
    family_ids: Optional[List[str]] = None,
) -> None:
    """
    Write a list of sample IDs to a file (for PLINK --keep/--remove).

    Args:
        sample_ids: List of sample IDs
        output_path: Output file path
        family_ids: Optional list of family IDs (defaults to sample IDs)
    """
    if family_ids is None:
        family_ids = sample_ids

    with open(output_path, "w") as f:
        for fid, iid in zip(family_ids, sample_ids):
            f.write(f"{fid}\t{iid}\n")


def read_snp_list(snp_file: Union[str, Path]) -> List[str]:
    """
    Read a list of SNP IDs from a file.

    Args:
        snp_file: Path to file containing SNP IDs (one per line)

    Returns:
        List of SNP IDs
    """
    snp_ids = []
    with open(snp_file, "r") as f:
        for line in f:
            snp_id = line.strip()
            if snp_id and not snp_id.startswith("#"):
                snp_ids.append(snp_id)
    return snp_ids


def read_sample_list(sample_file: Union[str, Path]) -> List[Tuple[str, str]]:
    """
    Read a list of sample IDs from a file (PLINK format: FID IID).

    Args:
        sample_file: Path to file containing sample IDs

    Returns:
        List of (family_id, individual_id) tuples
    """
    samples = []
    with open(sample_file, "r") as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 2:
                samples.append((parts[0], parts[1]))
            elif len(parts) == 1:
                samples.append((parts[0], parts[0]))
    return samples


if __name__ == "__main__":
    # Test basic functionality
    print("File Utils Test")
    print("=" * 50)

    # Test ensure_dir
    test_dir = ensure_dir("/tmp/test_ancestry_pipeline")
    print(f"Created directory: {test_dir}")

    # Test format_size
    print(f"1024 bytes = {format_size(1024)}")
    print(f"1048576 bytes = {format_size(1048576)}")
    print(f"1073741824 bytes = {format_size(1073741824)}")
