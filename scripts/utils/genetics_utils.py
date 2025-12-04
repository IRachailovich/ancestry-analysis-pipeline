#!/usr/bin/env python3
"""
Genetics-specific utilities for the ancestry analysis pipeline.

This module provides utilities for:
- SNP ID parsing and validation (rsID)
- Chromosome validation
- Allele strand flipping (A↔T, G↔C)
- Palindromic SNP detection
- Genotype encoding/decoding
- Minor allele frequency (MAF) calculation
- Hardy-Weinberg equilibrium testing

Example:
    >>> from scripts.utils.genetics_utils import is_valid_rsid, flip_allele
    >>> is_valid_rsid("rs12345")
    True
    >>> flip_allele("A")
    'T'
"""

import math
import re
from typing import Dict, List, Optional, Tuple, Union


# ==============================================================================
# Constants
# ==============================================================================

# Valid nucleotides
NUCLEOTIDES = {"A", "T", "G", "C"}

# Complement mapping for strand flipping
COMPLEMENT: Dict[str, str] = {
    "A": "T",
    "T": "A",
    "G": "C",
    "C": "G",
}

# Autosomal chromosome labels
AUTOSOMAL_CHROMOSOMES = {str(i) for i in range(1, 23)}
AUTOSOMAL_CHROMOSOMES_WITH_CHR = {f"chr{i}" for i in range(1, 23)}

# Sex and mitochondrial chromosomes
SEX_CHROMOSOMES = {"X", "Y", "chrX", "chrY", "23", "24"}
MITOCHONDRIAL = {"MT", "M", "chrM", "chrMT", "26"}

# All valid chromosomes
ALL_CHROMOSOMES = AUTOSOMAL_CHROMOSOMES | SEX_CHROMOSOMES | MITOCHONDRIAL

# rsID pattern
RSID_PATTERN = re.compile(r"^rs\d+$", re.IGNORECASE)


# ==============================================================================
# SNP ID Validation
# ==============================================================================

def is_valid_rsid(snp_id: str) -> bool:
    """
    Check if a SNP ID is a valid rsID format.

    Args:
        snp_id: SNP identifier to validate

    Returns:
        True if the ID matches rsID format (rs followed by digits)

    Example:
        >>> is_valid_rsid("rs12345")
        True
        >>> is_valid_rsid("chr1:12345")
        False
    """
    return bool(RSID_PATTERN.match(snp_id))


def parse_rsid(snp_id: str) -> Optional[int]:
    """
    Parse an rsID and return the numeric portion.

    Args:
        snp_id: SNP identifier

    Returns:
        The numeric ID, or None if not a valid rsID

    Example:
        >>> parse_rsid("rs12345")
        12345
    """
    if is_valid_rsid(snp_id):
        return int(snp_id[2:])
    return None


def normalize_snp_id(snp_id: str) -> str:
    """
    Normalize a SNP ID to lowercase rsID format.

    Args:
        snp_id: SNP identifier

    Returns:
        Normalized SNP ID

    Example:
        >>> normalize_snp_id("RS12345")
        'rs12345'
    """
    if is_valid_rsid(snp_id):
        return snp_id.lower()
    return snp_id


# ==============================================================================
# Chromosome Validation
# ==============================================================================

def is_autosomal(chromosome: Union[str, int]) -> bool:
    """
    Check if a chromosome is autosomal (1-22).

    Args:
        chromosome: Chromosome identifier (string or integer)

    Returns:
        True if the chromosome is autosomal

    Example:
        >>> is_autosomal("1")
        True
        >>> is_autosomal("X")
        False
        >>> is_autosomal(22)
        True
    """
    chrom_str = str(chromosome)

    # Handle 'chr' prefix
    if chrom_str.lower().startswith("chr"):
        chrom_str = chrom_str[3:]

    return chrom_str in AUTOSOMAL_CHROMOSOMES


def is_sex_chromosome(chromosome: Union[str, int]) -> bool:
    """
    Check if a chromosome is a sex chromosome (X or Y).

    Args:
        chromosome: Chromosome identifier

    Returns:
        True if the chromosome is X or Y
    """
    return str(chromosome) in SEX_CHROMOSOMES


def is_mitochondrial(chromosome: Union[str, int]) -> bool:
    """
    Check if a chromosome is mitochondrial.

    Args:
        chromosome: Chromosome identifier

    Returns:
        True if the chromosome is mitochondrial
    """
    return str(chromosome) in MITOCHONDRIAL


def normalize_chromosome(chromosome: Union[str, int]) -> str:
    """
    Normalize chromosome notation to simple format (1-22, X, Y, MT).

    Args:
        chromosome: Chromosome identifier

    Returns:
        Normalized chromosome string

    Example:
        >>> normalize_chromosome("chr1")
        '1'
        >>> normalize_chromosome("chrX")
        'X'
        >>> normalize_chromosome(22)
        '22'
    """
    chrom_str = str(chromosome)

    # Remove 'chr' prefix
    if chrom_str.lower().startswith("chr"):
        chrom_str = chrom_str[3:]

    # Normalize sex chromosomes
    if chrom_str.upper() in {"X", "23"}:
        return "X"
    if chrom_str.upper() in {"Y", "24"}:
        return "Y"
    if chrom_str.upper() in {"M", "MT", "26"}:
        return "MT"

    return chrom_str


def chromosome_sort_key(chromosome: Union[str, int]) -> Tuple[int, str]:
    """
    Get a sort key for chromosome ordering.

    Orders: 1-22, X, Y, MT

    Args:
        chromosome: Chromosome identifier

    Returns:
        Tuple for sorting (numeric order, string fallback)
    """
    chrom = normalize_chromosome(chromosome)

    # Try to convert to integer
    try:
        return (int(chrom), "")
    except ValueError:
        pass

    # Handle sex and mitochondrial chromosomes
    special_order = {"X": 23, "Y": 24, "MT": 25}
    if chrom.upper() in special_order:
        return (special_order[chrom.upper()], "")

    # Fallback
    return (100, chrom)


# ==============================================================================
# Allele Manipulation
# ==============================================================================

def is_valid_allele(allele: str) -> bool:
    """
    Check if an allele is a valid single nucleotide.

    Args:
        allele: Allele string

    Returns:
        True if the allele is A, T, G, or C
    """
    return allele.upper() in NUCLEOTIDES


def flip_allele(allele: str) -> str:
    """
    Get the complement of an allele (strand flip).

    Args:
        allele: Single nucleotide (A, T, G, or C)

    Returns:
        Complementary nucleotide

    Raises:
        ValueError: If the allele is not a valid nucleotide

    Example:
        >>> flip_allele("A")
        'T'
        >>> flip_allele("G")
        'C'
    """
    allele_upper = allele.upper()
    if allele_upper not in COMPLEMENT:
        raise ValueError(f"Invalid allele: {allele}. Must be A, T, G, or C.")
    return COMPLEMENT[allele_upper]


def flip_genotype(genotype: str) -> str:
    """
    Flip both alleles in a genotype string.

    Args:
        genotype: Two-character genotype (e.g., "AT", "GC")

    Returns:
        Flipped genotype

    Example:
        >>> flip_genotype("AT")
        'TA'
        >>> flip_genotype("GC")
        'CG'
    """
    if len(genotype) != 2:
        raise ValueError(f"Genotype must be exactly 2 characters: {genotype}")

    return flip_allele(genotype[0]) + flip_allele(genotype[1])


def is_palindromic(allele1: str, allele2: str) -> bool:
    """
    Check if a SNP is palindromic (A/T or G/C).

    Palindromic SNPs have alleles that are complements of each other,
    making strand orientation ambiguous.

    Args:
        allele1: First allele
        allele2: Second allele

    Returns:
        True if the SNP is palindromic

    Example:
        >>> is_palindromic("A", "T")
        True
        >>> is_palindromic("G", "C")
        True
        >>> is_palindromic("A", "G")
        False
    """
    a1, a2 = allele1.upper(), allele2.upper()

    # Check both orderings
    palindromic_pairs = {("A", "T"), ("T", "A"), ("G", "C"), ("C", "G")}
    return (a1, a2) in palindromic_pairs


def are_alleles_compatible(
    alleles1: Tuple[str, str],
    alleles2: Tuple[str, str],
) -> Tuple[bool, bool]:
    """
    Check if two sets of alleles are compatible (possibly on opposite strands).

    Args:
        alleles1: First (ref, alt) allele pair
        alleles2: Second (ref, alt) allele pair

    Returns:
        Tuple of (compatible, needs_flip)
        - compatible: True if alleles match (directly or with strand flip)
        - needs_flip: True if strand flip is needed for matching

    Example:
        >>> are_alleles_compatible(("A", "G"), ("A", "G"))
        (True, False)
        >>> are_alleles_compatible(("A", "G"), ("T", "C"))
        (True, True)
    """
    a1_ref, a1_alt = alleles1[0].upper(), alleles1[1].upper()
    a2_ref, a2_alt = alleles2[0].upper(), alleles2[1].upper()

    # Direct match
    if (a1_ref == a2_ref and a1_alt == a2_alt) or (a1_ref == a2_alt and a1_alt == a2_ref):
        return (True, False)

    # Strand flip match
    try:
        a1_ref_flip = flip_allele(a1_ref)
        a1_alt_flip = flip_allele(a1_alt)

        if (a1_ref_flip == a2_ref and a1_alt_flip == a2_alt) or \
           (a1_ref_flip == a2_alt and a1_alt_flip == a2_ref):
            return (True, True)
    except ValueError:
        pass

    return (False, False)


# ==============================================================================
# Genotype Encoding/Decoding
# ==============================================================================

def encode_genotype_012(allele1: str, allele2: str, ref_allele: str) -> int:
    """
    Encode a genotype as 0/1/2 based on reference allele count.

    Args:
        allele1: First allele
        allele2: Second allele
        ref_allele: Reference allele

    Returns:
        0 (homozygous alt), 1 (heterozygous), or 2 (homozygous ref)
        Returns -1 for missing genotypes

    Example:
        >>> encode_genotype_012("A", "A", "A")
        2
        >>> encode_genotype_012("A", "G", "A")
        1
        >>> encode_genotype_012("G", "G", "A")
        0
    """
    # Handle missing
    if allele1 in {"0", "-", "N", "."} or allele2 in {"0", "-", "N", "."}:
        return -1

    a1, a2, ref = allele1.upper(), allele2.upper(), ref_allele.upper()

    ref_count = (1 if a1 == ref else 0) + (1 if a2 == ref else 0)
    return ref_count


def decode_012_genotype(
    code: int,
    ref_allele: str,
    alt_allele: str,
) -> Tuple[str, str]:
    """
    Decode a 0/1/2 genotype code to alleles.

    Args:
        code: Genotype code (0, 1, or 2)
        ref_allele: Reference allele
        alt_allele: Alternate allele

    Returns:
        Tuple of (allele1, allele2)

    Raises:
        ValueError: If code is not 0, 1, or 2
    """
    if code == 2:
        return (ref_allele, ref_allele)
    elif code == 1:
        return (ref_allele, alt_allele)
    elif code == 0:
        return (alt_allele, alt_allele)
    else:
        raise ValueError(f"Invalid genotype code: {code}. Must be 0, 1, or 2.")


def genotype_to_eigenstrat(allele1: str, allele2: str, ref_allele: str) -> str:
    """
    Convert genotype to EIGENSTRAT format (0/1/2/9).

    Args:
        allele1: First allele
        allele2: Second allele
        ref_allele: Reference allele

    Returns:
        EIGENSTRAT genotype code as string
    """
    code = encode_genotype_012(allele1, allele2, ref_allele)
    return "9" if code == -1 else str(code)


# ==============================================================================
# Minor Allele Frequency
# ==============================================================================

def calculate_maf(
    genotype_counts: Tuple[int, int, int],
) -> float:
    """
    Calculate minor allele frequency from genotype counts.

    Args:
        genotype_counts: Tuple of (count_0, count_1, count_2) where:
            - count_0: Number of homozygous alternate individuals
            - count_1: Number of heterozygous individuals
            - count_2: Number of homozygous reference individuals

    Returns:
        Minor allele frequency (0.0 to 0.5)

    Example:
        >>> calculate_maf((10, 40, 50))  # 100 individuals
        0.3
    """
    n0, n1, n2 = genotype_counts
    total = n0 + n1 + n2

    if total == 0:
        return 0.0

    # Frequency of alternate allele
    alt_freq = (2 * n0 + n1) / (2 * total)

    # MAF is the lesser of ref and alt frequency
    return min(alt_freq, 1 - alt_freq)


def calculate_maf_from_genotypes(genotypes: List[int]) -> float:
    """
    Calculate minor allele frequency from a list of genotypes.

    Args:
        genotypes: List of genotype codes (0, 1, 2, or -1 for missing)

    Returns:
        Minor allele frequency
    """
    # Filter out missing genotypes
    valid_genotypes = [g for g in genotypes if g >= 0]

    if not valid_genotypes:
        return 0.0

    n0 = sum(1 for g in valid_genotypes if g == 0)
    n1 = sum(1 for g in valid_genotypes if g == 1)
    n2 = sum(1 for g in valid_genotypes if g == 2)

    return calculate_maf((n0, n1, n2))


# ==============================================================================
# Hardy-Weinberg Equilibrium
# ==============================================================================

def hardy_weinberg_test(
    genotype_counts: Tuple[int, int, int],
) -> Tuple[float, float]:
    """
    Perform Hardy-Weinberg equilibrium exact test.

    Uses the chi-squared approximation for large sample sizes.

    Args:
        genotype_counts: Tuple of (count_0, count_1, count_2) where:
            - count_0: Number of homozygous alternate (aa)
            - count_1: Number of heterozygous (Aa)
            - count_2: Number of homozygous reference (AA)

    Returns:
        Tuple of (chi_squared, p_value)

    Example:
        >>> chi2, pval = hardy_weinberg_test((10, 50, 40))
        >>> print(f"Chi-squared: {chi2:.2f}, p-value: {pval:.4f}")
    """
    n0, n1, n2 = genotype_counts
    n = n0 + n1 + n2

    if n == 0:
        return (0.0, 1.0)

    # Observed allele frequencies
    p = (2 * n2 + n1) / (2 * n)  # Frequency of reference allele
    q = 1 - p  # Frequency of alternate allele

    # Expected counts under HWE
    e2 = n * p * p  # Expected AA
    e1 = n * 2 * p * q  # Expected Aa
    e0 = n * q * q  # Expected aa

    # Chi-squared statistic
    chi2 = 0.0
    for obs, exp in [(n2, e2), (n1, e1), (n0, e0)]:
        if exp > 0:
            chi2 += (obs - exp) ** 2 / exp

    # P-value from chi-squared distribution (df=1)
    # Using approximation: P(X > chi2) for chi-squared with 1 df
    pval = chi_squared_pvalue(chi2, df=1)

    return (chi2, pval)


def chi_squared_pvalue(chi2: float, df: int = 1) -> float:
    """
    Calculate p-value for chi-squared test.

    Uses the regularized incomplete gamma function approximation.

    Args:
        chi2: Chi-squared statistic
        df: Degrees of freedom

    Returns:
        P-value
    """
    if chi2 <= 0:
        return 1.0

    # Using the incomplete gamma function approximation
    # P-value = 1 - regularized_gamma_lower(df/2, chi2/2)
    # For df=1, this simplifies to 2 * (1 - erf(sqrt(chi2/2)))

    if df == 1:
        # Special case for df=1 using error function approximation
        x = math.sqrt(chi2 / 2)
        return 2 * (1 - _erf_approx(x))
    else:
        # General case using Gamma function approximation
        return _incomplete_gamma_upper(df / 2, chi2 / 2)


def _erf_approx(x: float) -> float:
    """
    Approximate error function using Abramowitz and Stegun method.
    """
    # Constants
    a1 = 0.254829592
    a2 = -0.284496736
    a3 = 1.421413741
    a4 = -1.453152027
    a5 = 1.061405429
    p = 0.3275911

    sign = 1 if x >= 0 else -1
    x = abs(x)

    t = 1.0 / (1.0 + p * x)
    y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * math.exp(-x * x)

    return sign * y


def _incomplete_gamma_upper(a: float, x: float) -> float:
    """
    Approximate upper incomplete gamma function Q(a, x).

    Uses series expansion for small x and continued fraction for large x.
    """
    if x < 0 or a <= 0:
        return 1.0

    if x < a + 1:
        # Use series expansion
        return 1.0 - _gamma_series(a, x)
    else:
        # Use continued fraction
        return _gamma_cf(a, x)


def _gamma_series(a: float, x: float, max_iter: int = 100, eps: float = 1e-10) -> float:
    """Series expansion for lower incomplete gamma function."""
    if x == 0:
        return 0.0

    ap = a
    sum_val = 1.0 / a
    delta = sum_val

    for _ in range(max_iter):
        ap += 1
        delta *= x / ap
        sum_val += delta
        if abs(delta) < abs(sum_val) * eps:
            break

    return sum_val * math.exp(-x + a * math.log(x) - math.lgamma(a))


def _gamma_cf(a: float, x: float, max_iter: int = 100, eps: float = 1e-10) -> float:
    """Continued fraction for upper incomplete gamma function."""
    b = x + 1 - a
    c = 1.0 / 1e-30
    d = 1.0 / b
    h = d

    for i in range(1, max_iter + 1):
        an = -i * (i - a)
        b += 2
        d = an * d + b
        if abs(d) < 1e-30:
            d = 1e-30
        c = b + an / c
        if abs(c) < 1e-30:
            c = 1e-30
        d = 1.0 / d
        delta = d * c
        h *= delta
        if abs(delta - 1) < eps:
            break

    return h * math.exp(-x + a * math.log(x) - math.lgamma(a))


# ==============================================================================
# 23andMe Format Utilities
# ==============================================================================

def detect_23andme_version(header_lines: List[str]) -> str:
    """
    Detect 23andMe file format version from header comments.

    23andMe has released several file format versions over the years:
    - v3: Uses NCBI Build 36 (hg18), older chip versions
    - v4: Uses GRCh37/hg19 (Build 37), introduced around 2013-2017
    - v5: Uses GRCh38/hg38 (Build 38), current format since ~2017

    The function detects the version by looking for build information in
    the header comments. If no explicit build info is found but "23andme"
    is mentioned, it defaults to v5 as the most recent format.

    Args:
        header_lines: List of header lines (starting with #)

    Returns:
        Version string ("v3", "v4", "v5", or "unknown")

    Example:
        >>> detect_23andme_version(["# build 37"])
        'v4'
        >>> detect_23andme_version(["# This data file generated by 23andMe..."])
        'v5'
    """
    header_text = "\n".join(header_lines).lower()

    if "build 37" in header_text or "grch37" in header_text:
        return "v4"
    elif "build 38" in header_text or "grch38" in header_text:
        return "v5"
    elif "build 36" in header_text or "ncbi36" in header_text:
        return "v3"
    elif "23andme" in header_text:
        # Assume v5 for newer files without explicit build info
        return "v5"

    return "unknown"


def parse_23andme_genotype(genotype: str) -> Tuple[str, str]:
    """
    Parse a 23andMe genotype string.

    Args:
        genotype: Genotype string (e.g., "AA", "AG", "--", "DD", "II")

    Returns:
        Tuple of (allele1, allele2)

    Example:
        >>> parse_23andme_genotype("AG")
        ('A', 'G')
        >>> parse_23andme_genotype("--")
        ('0', '0')
    """
    genotype = genotype.strip().upper()

    # Handle missing genotypes
    if genotype in {"--", "NN", "00", ""}:
        return ("0", "0")

    # Handle insertions/deletions
    if genotype in {"DD", "DI", "ID", "II"}:
        return ("0", "0")

    # Normal two-allele genotype
    if len(genotype) == 2:
        return (genotype[0], genotype[1])

    # Single character (some formats use this for hemizygous calls on X/Y)
    if len(genotype) == 1:
        return (genotype, genotype)

    return ("0", "0")


def is_valid_23andme_snp(
    chromosome: str,
    position: int,
    genotype: str,
    filter_autosomal: bool = True,
) -> bool:
    """
    Check if a 23andMe SNP line is valid for analysis.

    Args:
        chromosome: Chromosome string
        position: Base pair position
        genotype: Genotype string
        filter_autosomal: If True, only accept autosomal SNPs

    Returns:
        True if the SNP is valid
    """
    # Check genotype
    a1, a2 = parse_23andme_genotype(genotype)
    if a1 == "0" or a2 == "0":
        return False

    # Check position
    if position <= 0:
        return False

    # Check chromosome
    if filter_autosomal and not is_autosomal(chromosome):
        return False

    return True


if __name__ == "__main__":
    # Test basic functionality
    print("Genetics Utils Test")
    print("=" * 50)

    # Test rsID validation
    print(f"is_valid_rsid('rs12345'): {is_valid_rsid('rs12345')}")
    print(f"is_valid_rsid('chr1:12345'): {is_valid_rsid('chr1:12345')}")

    # Test chromosome validation
    print(f"is_autosomal('1'): {is_autosomal('1')}")
    print(f"is_autosomal('X'): {is_autosomal('X')}")

    # Test allele flipping
    print(f"flip_allele('A'): {flip_allele('A')}")
    print(f"flip_allele('G'): {flip_allele('G')}")

    # Test palindromic detection
    print(f"is_palindromic('A', 'T'): {is_palindromic('A', 'T')}")
    print(f"is_palindromic('A', 'G'): {is_palindromic('A', 'G')}")

    # Test MAF calculation
    print(f"calculate_maf((10, 40, 50)): {calculate_maf((10, 40, 50))}")

    # Test HWE
    chi2, pval = hardy_weinberg_test((10, 50, 40))
    print(f"HWE test (10, 50, 40): chi2={chi2:.2f}, p={pval:.4f}")
