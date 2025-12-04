#!/usr/bin/env python3
"""
Configuration settings for the ancestry analysis pipeline.

This module provides centralized configuration for:
- Directory paths (data, results, logs)
- External tool paths (PLINK, ADMIXTURE, EIGENSOFT, ADMIXTOOLS)
- Reference population definitions and colors
- Analysis parameters (MAF, LD, missingness thresholds)

All paths are configurable via environment variables.

Example:
    >>> from scripts.utils.config import PROJECT_ROOT, DATA_DIR
    >>> from scripts.utils.config import get_tool_path
    >>> plink_path = get_tool_path("plink")
"""

import os
import shutil
from pathlib import Path
from typing import Dict, List, Optional


# ==============================================================================
# Project Directory Structure
# ==============================================================================

def _get_project_root() -> Path:
    """Determine the project root directory."""
    env_root = os.environ.get("ANCESTRY_PROJECT_ROOT")
    if env_root:
        return Path(env_root).resolve()
    # Default: assume this file is in scripts/utils/
    return Path(__file__).resolve().parent.parent.parent


PROJECT_ROOT: Path = _get_project_root()
"""Project root directory."""

DATA_DIR: Path = Path(os.environ.get("ANCESTRY_DATA_DIR", PROJECT_ROOT / "data"))
"""Base data directory."""

RAW_DATA_DIR: Path = DATA_DIR / "raw"
"""Raw input data directory."""

PROCESSED_DATA_DIR: Path = DATA_DIR / "processed"
"""Processed data directory."""

RESULTS_DIR: Path = Path(os.environ.get("ANCESTRY_RESULTS_DIR", PROJECT_ROOT / "results"))
"""Analysis results directory."""

LOGS_DIR: Path = Path(os.environ.get("ANCESTRY_LOGS_DIR", PROJECT_ROOT / "logs"))
"""Log files directory."""


# ==============================================================================
# External Tool Paths
# ==============================================================================

# Default tool installation directory
DEFAULT_TOOLS_DIR: Path = Path(os.environ.get("ANCESTRY_TOOLS_DIR", Path.home() / "tools"))

# Tool-specific paths (configurable via environment variables)
TOOL_PATHS: Dict[str, Optional[str]] = {
    "plink": os.environ.get("PLINK_PATH"),
    "plink2": os.environ.get("PLINK2_PATH"),
    "admixture": os.environ.get("ADMIXTURE_PATH"),
    "smartpca": os.environ.get("SMARTPCA_PATH"),
    "convertf": os.environ.get("CONVERTF_PATH"),
    "qpadm": os.environ.get("QPADM_PATH"),
    "qpdstat": os.environ.get("QPDSTAT_PATH"),
    "qpgraph": os.environ.get("QPGRAPH_PATH"),
    "rfmix": os.environ.get("RFMIX_PATH"),
    "treemix": os.environ.get("TREEMIX_PATH"),
}


def get_tool_path(tool_name: str) -> Optional[Path]:
    """
    Get the path to an external tool.

    First checks environment variables, then looks in the default tools
    directory, and finally checks if the tool is in PATH.

    Args:
        tool_name: Name of the tool (e.g., "plink", "admixture", "smartpca")

    Returns:
        Path to the tool executable, or None if not found

    Raises:
        ValueError: If tool_name is not a recognized tool

    Example:
        >>> plink = get_tool_path("plink")
        >>> if plink:
        ...     print(f"PLINK found at: {plink}")
    """
    tool_name_lower = tool_name.lower()

    if tool_name_lower not in TOOL_PATHS:
        raise ValueError(
            f"Unknown tool: {tool_name}. "
            f"Known tools: {', '.join(TOOL_PATHS.keys())}"
        )

    # Check environment variable first
    env_path = TOOL_PATHS.get(tool_name_lower)
    if env_path:
        path = Path(env_path)
        if path.exists() and path.is_file():
            return path

    # Check default tools directory
    tool_locations = [
        DEFAULT_TOOLS_DIR / tool_name_lower / tool_name_lower,
        DEFAULT_TOOLS_DIR / tool_name_lower / f"{tool_name_lower}.exe",
        DEFAULT_TOOLS_DIR / tool_name_lower,
        DEFAULT_TOOLS_DIR / f"{tool_name_lower}.exe",
    ]

    for loc in tool_locations:
        if loc.exists() and loc.is_file():
            return loc

    # Check if tool is in PATH
    which_result = shutil.which(tool_name_lower)
    if which_result:
        return Path(which_result)

    return None


def check_tool_available(tool_name: str) -> bool:
    """
    Check if an external tool is available.

    Args:
        tool_name: Name of the tool to check

    Returns:
        True if the tool is available, False otherwise
    """
    return get_tool_path(tool_name) is not None


def get_required_tools() -> List[str]:
    """
    Get list of tools required for the full pipeline.

    Returns:
        List of tool names
    """
    return ["plink", "admixture", "smartpca", "convertf"]


def check_all_tools() -> Dict[str, bool]:
    """
    Check availability of all registered tools.

    Returns:
        Dictionary mapping tool names to availability status
    """
    return {tool: check_tool_available(tool) for tool in TOOL_PATHS.keys()}


# ==============================================================================
# Reference Population Definitions
# ==============================================================================

# Population colors for plotting
POPULATION_COLORS: Dict[str, str] = {
    # 1000 Genomes superpopulations
    "AFR": "#E41A1C",  # African - Red
    "EUR": "#377EB8",  # European - Blue
    "EAS": "#4DAF4A",  # East Asian - Green
    "SAS": "#984EA3",  # South Asian - Purple
    "AMR": "#FF7F00",  # American - Orange

    # 1000 Genomes populations
    "YRI": "#E41A1C",
    "LWK": "#B2182B",
    "GWD": "#D6604D",
    "MSL": "#F4A582",
    "ESN": "#FDDBC7",
    "ASW": "#D1E5F0",
    "ACB": "#92C5DE",
    "CEU": "#377EB8",
    "TSI": "#4393C3",
    "FIN": "#2166AC",
    "GBR": "#053061",
    "IBS": "#67A9CF",
    "CHB": "#4DAF4A",
    "JPT": "#1B7837",
    "CHS": "#5AAE61",
    "CDX": "#A6D96A",
    "KHV": "#D9F0D3",
    "GIH": "#984EA3",
    "PJL": "#762A83",
    "BEB": "#9970AB",
    "STU": "#C2A5CF",
    "ITU": "#E7D4E8",
    "MXL": "#FF7F00",
    "PUR": "#E66101",
    "CLM": "#FDB863",
    "PEL": "#FEE0B6",

    # Ancient populations
    "Yamnaya": "#8B4513",  # Saddle Brown
    "Yamnaya_Samara": "#8B4513",
    "Anatolia_N": "#228B22",  # Forest Green
    "WHG": "#4169E1",  # Royal Blue
    "EHG": "#DC143C",  # Crimson
    "CHG": "#9932CC",  # Dark Orchid
    "Steppe": "#A0522D",  # Sienna
    "Farmer": "#32CD32",  # Lime Green
    "HG": "#1E90FF",  # Dodger Blue

    # Sample
    "Sample": "#000000",  # Black
    "TARGET": "#000000",
}

# Superpopulation mappings
SUPERPOPULATION_MAP: Dict[str, str] = {
    "YRI": "AFR", "LWK": "AFR", "GWD": "AFR", "MSL": "AFR", "ESN": "AFR",
    "ASW": "AFR", "ACB": "AFR",
    "CEU": "EUR", "TSI": "EUR", "FIN": "EUR", "GBR": "EUR", "IBS": "EUR",
    "CHB": "EAS", "JPT": "EAS", "CHS": "EAS", "CDX": "EAS", "KHV": "EAS",
    "GIH": "SAS", "PJL": "SAS", "BEB": "SAS", "STU": "SAS", "ITU": "SAS",
    "MXL": "AMR", "PUR": "AMR", "CLM": "AMR", "PEL": "AMR",
}


def get_population_color(population: str) -> str:
    """
    Get the plotting color for a population.

    Args:
        population: Population label

    Returns:
        Hex color code
    """
    return POPULATION_COLORS.get(population, "#808080")  # Default gray


# ==============================================================================
# Analysis Parameters
# ==============================================================================

# Quality Control Defaults
QC_DEFAULTS: Dict[str, float] = {
    "snp_missingness": 0.05,  # Maximum SNP missingness rate
    "ind_missingness": 0.10,  # Maximum individual missingness rate
    "maf": 0.01,  # Minimum minor allele frequency
    "hwe_p": 1e-6,  # Hardy-Weinberg equilibrium p-value threshold
}

# LD Pruning Defaults
LD_DEFAULTS: Dict[str, float] = {
    "window_kb": 50,  # Window size in kb
    "step": 5,  # Step size (SNPs)
    "r2_threshold": 0.2,  # r² threshold
}

# Relatedness Defaults
RELATEDNESS_DEFAULTS: Dict[str, float] = {
    "pi_hat": 0.2,  # PI_HAT threshold for relatedness
}

# PCA Defaults
PCA_DEFAULTS: Dict[str, int] = {
    "n_pcs": 10,  # Number of principal components
    "n_outlier_iter": 5,  # Number of outlier removal iterations
}

# ADMIXTURE Defaults
ADMIXTURE_DEFAULTS: Dict[str, int] = {
    "min_k": 2,  # Minimum K value
    "max_k": 15,  # Maximum K value
    "cv_folds": 5,  # Cross-validation folds
    "threads": 4,  # Number of threads
}


# ==============================================================================
# Ancient DNA Source Populations
# ==============================================================================

# Common ancient DNA source populations for qpAdm
ANCIENT_SOURCES: Dict[str, List[str]] = {
    # Western Hunter-Gatherers
    "WHG": ["Loschbour", "Bichon", "La_Brana1"],

    # Eastern Hunter-Gatherers
    "EHG": ["Karelia_HG", "Samara_HG"],

    # Caucasus Hunter-Gatherers
    "CHG": ["Kotias", "Satsurblia"],

    # Anatolian Neolithic
    "Anatolia_N": ["Barcin_N", "Mentese_N"],

    # Yamnaya Steppe pastoralists
    "Yamnaya": ["Yamnaya_Samara", "Yamnaya_Kalmykia"],

    # Standard outgroups for qpAdm
    "Outgroups": [
        "Mbuti.DG",
        "Papuan.DG",
        "Onge.DG",
        "Han.DG",
        "Karitiana.DG",
    ],
}

# Standard qpAdm model configurations
QPADM_MODELS: Dict[str, Dict[str, List[str]]] = {
    "european_2way": {
        "sources": ["Yamnaya_Samara", "Anatolia_N"],
        "outgroups": ["Mbuti.DG", "Papuan.DG", "Onge.DG", "Han.DG"],
    },
    "european_3way": {
        "sources": ["Yamnaya_Samara", "Anatolia_N", "WHG"],
        "outgroups": ["Mbuti.DG", "Papuan.DG", "Onge.DG", "Han.DG", "Karitiana.DG"],
    },
}


# ==============================================================================
# Reference Dataset Paths
# ==============================================================================

def get_reference_path(dataset: str) -> Path:
    """
    Get the path to a reference dataset.

    Args:
        dataset: Dataset name ("1000g", "aadr", "sgdp")

    Returns:
        Path to the dataset directory

    Raises:
        ValueError: If dataset name is not recognized
    """
    dataset_map = {
        "1000g": RAW_DATA_DIR / "1000genomes",
        "1000genomes": RAW_DATA_DIR / "1000genomes",
        "aadr": RAW_DATA_DIR / "aadr",
        "sgdp": RAW_DATA_DIR / "sgdp",
    }

    dataset_lower = dataset.lower()
    if dataset_lower not in dataset_map:
        raise ValueError(
            f"Unknown dataset: {dataset}. "
            f"Known datasets: 1000g, aadr, sgdp"
        )

    return dataset_map[dataset_lower]


# ==============================================================================
# Logging Configuration
# ==============================================================================

LOG_FORMAT: str = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
LOG_DATE_FORMAT: str = "%Y-%m-%d %H:%M:%S"


def setup_logging(
    name: str,
    log_file: Optional[Path] = None,
    level: int = 20,  # logging.INFO
) -> None:
    """
    Configure logging for a script.

    Args:
        name: Logger name
        log_file: Path to log file (optional)
        level: Logging level (default: INFO)
    """
    import logging

    logger = logging.getLogger(name)
    logger.setLevel(level)

    # Console handler
    console_handler = logging.StreamHandler()
    console_handler.setLevel(level)
    console_formatter = logging.Formatter("%(levelname)s: %(message)s")
    console_handler.setFormatter(console_formatter)
    logger.addHandler(console_handler)

    # File handler (if specified)
    if log_file:
        log_file.parent.mkdir(parents=True, exist_ok=True)
        file_handler = logging.FileHandler(log_file)
        file_handler.setLevel(logging.DEBUG)
        file_formatter = logging.Formatter(LOG_FORMAT, datefmt=LOG_DATE_FORMAT)
        file_handler.setFormatter(file_formatter)
        logger.addHandler(file_handler)


if __name__ == "__main__":
    # Print configuration when run directly
    print("Ancestry Analysis Pipeline Configuration")
    print("=" * 50)
    print(f"Project Root: {PROJECT_ROOT}")
    print(f"Data Directory: {DATA_DIR}")
    print(f"Results Directory: {RESULTS_DIR}")
    print(f"Logs Directory: {LOGS_DIR}")
    print()
    print("Tool Availability:")
    for tool, available in check_all_tools().items():
        status = "✓" if available else "✗"
        path = get_tool_path(tool) or "Not found"
        print(f"  {status} {tool}: {path}")
