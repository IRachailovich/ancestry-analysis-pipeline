"""
Utility modules for the ancestry analysis pipeline.

This package provides shared utilities for:
- Configuration management (config.py)
- File handling for genetic data formats (file_utils.py)
- Genetics-specific operations (genetics_utils.py)
"""

from .config import (
    PROJECT_ROOT,
    DATA_DIR,
    RESULTS_DIR,
    LOGS_DIR,
    get_tool_path,
)
from .file_utils import (
    check_plink_files,
    check_eigenstrat_files,
    ensure_dir,
)
from .genetics_utils import (
    is_valid_rsid,
    is_autosomal,
    flip_allele,
    is_palindromic,
)

__all__ = [
    "PROJECT_ROOT",
    "DATA_DIR",
    "RESULTS_DIR",
    "LOGS_DIR",
    "get_tool_path",
    "check_plink_files",
    "check_eigenstrat_files",
    "ensure_dir",
    "is_valid_rsid",
    "is_autosomal",
    "flip_allele",
    "is_palindromic",
]
