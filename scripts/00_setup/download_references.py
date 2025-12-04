#!/usr/bin/env python3
"""
Download reference datasets for ancestry analysis.

Downloads:
- 1000 Genomes Project Phase 3
- Reich Lab Ancient DNA Resource (AADR)
- Simons Genome Diversity Project (SGDP)

Usage:
    python download_references.py --dataset all
    python download_references.py --dataset 1000g
    python download_references.py --dataset aadr --force
    python download_references.py --dataset sgdp --output-dir /path/to/data/

Estimated download sizes:
    - 1000 Genomes: ~15 GB (VCF files for all chromosomes)
    - AADR: ~2 GB (v54.1 1240k dataset)
    - SGDP: ~500 MB (EIGENSTRAT format)
"""

import argparse
import ftplib
import logging
import os
import sys
import time
from pathlib import Path
from typing import Optional, List, Dict, Any
from urllib.parse import urlparse

import requests
from tqdm import tqdm


# Dataset URLs and file information
DATASET_INFO = {
    "1000g": {
        "description": "1000 Genomes Project Phase 3",
        "base_url": "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/",
        "files": [
            f"ALL.chr{i}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
            for i in range(1, 23)
        ],
        "panel_file": "integrated_call_samples_v3.20130502.ALL.panel",
        "output_subdir": "1000genomes",
        "estimated_size_gb": 15,
    },
    "aadr": {
        "description": "Reich Lab Ancient DNA Resource (AADR)",
        "base_url": "https://reichdata.hms.harvard.edu/pub/datasets/amh_repo/curated_releases/V54/V54.1/SHARE/public.dir/",
        "files": [
            "v54.1_1240K_public.geno",
            "v54.1_1240K_public.snp",
            "v54.1_1240K_public.ind",
        ],
        "output_subdir": "aadr",
        "estimated_size_gb": 2,
    },
    "sgdp": {
        "description": "Simons Genome Diversity Project",
        "base_url": "https://reichdata.hms.harvard.edu/pub/datasets/sgdp/",
        "files": [
            "cteam_extended.v4.geno",
            "cteam_extended.v4.snp",
            "cteam_extended.v4.ind",
        ],
        "output_subdir": "sgdp",
        "estimated_size_gb": 0.5,
    },
}


def setup_logging(log_dir: Path) -> logging.Logger:
    """Set up logging to console and file."""
    log_dir.mkdir(parents=True, exist_ok=True)
    log_file = log_dir / "download.log"

    # Create logger
    logger = logging.getLogger("download_references")
    logger.setLevel(logging.DEBUG)

    # Create formatters
    console_formatter = logging.Formatter("%(levelname)s: %(message)s")
    file_formatter = logging.Formatter(
        "%(asctime)s - %(levelname)s - %(message)s", datefmt="%Y-%m-%d %H:%M:%S"
    )

    # Console handler
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(logging.INFO)
    console_handler.setFormatter(console_formatter)

    # File handler
    file_handler = logging.FileHandler(log_file)
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(file_formatter)

    # Add handlers
    logger.addHandler(console_handler)
    logger.addHandler(file_handler)

    return logger


def get_file_size_ftp(host: str, path: str) -> Optional[int]:
    """Get file size from FTP server."""
    try:
        with ftplib.FTP(host) as ftp:
            ftp.login()
            return ftp.size(path)
    except Exception:
        return None


def get_file_size_http(url: str) -> Optional[int]:
    """Get file size from HTTP server."""
    try:
        response = requests.head(url, allow_redirects=True, timeout=30)
        if response.status_code == 200:
            content_length = response.headers.get("Content-Length")
            if content_length:
                return int(content_length)
    except Exception:
        pass
    return None


def format_size(size_bytes: int) -> str:
    """Format file size in human-readable format."""
    for unit in ["B", "KB", "MB", "GB", "TB"]:
        if size_bytes < 1024.0:
            return f"{size_bytes:.2f} {unit}"
        size_bytes /= 1024.0
    return f"{size_bytes:.2f} PB"


def download_file_http(
    url: str,
    output_path: Path,
    force: bool = False,
    max_retries: int = 3,
    logger: Optional[logging.Logger] = None,
) -> bool:
    """
    Download a file via HTTP with progress bar and error handling.

    Args:
        url: URL to download from
        output_path: Path to save the file
        force: Force re-download even if file exists
        max_retries: Number of retry attempts
        logger: Logger instance

    Returns:
        True if download successful, False otherwise
    """
    if logger is None:
        logger = logging.getLogger("download_references")

    # Check if file already exists
    if output_path.exists() and not force:
        logger.info(f"File already exists, skipping: {output_path.name}")
        return True

    # Create parent directories
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Get file size for progress bar
    expected_size = get_file_size_http(url)

    # Check for partial download (resume support)
    resume_pos = 0
    temp_path = output_path.with_suffix(output_path.suffix + ".partial")
    if temp_path.exists():
        resume_pos = temp_path.stat().st_size
        logger.info(f"Resuming download from {format_size(resume_pos)}")

    for attempt in range(1, max_retries + 1):
        try:
            headers = {}
            if resume_pos > 0:
                headers["Range"] = f"bytes={resume_pos}-"

            response = requests.get(
                url, headers=headers, stream=True, timeout=60, allow_redirects=True
            )

            # Handle range request response
            if resume_pos > 0 and response.status_code == 206:
                mode = "ab"
                total_size = resume_pos + int(
                    response.headers.get("Content-Length", 0)
                )
            elif response.status_code == 200:
                mode = "wb"
                resume_pos = 0
                total_size = int(response.headers.get("Content-Length", 0))
            else:
                logger.error(
                    f"HTTP error {response.status_code} for {url}"
                )
                continue

            # Download with progress bar
            desc = output_path.name[:40]
            with open(temp_path, mode) as f:
                with tqdm(
                    total=total_size,
                    initial=resume_pos,
                    unit="B",
                    unit_scale=True,
                    unit_divisor=1024,
                    desc=desc,
                    leave=True,
                ) as pbar:
                    for chunk in response.iter_content(chunk_size=8192):
                        if chunk:
                            f.write(chunk)
                            pbar.update(len(chunk))

            # Verify file size
            downloaded_size = temp_path.stat().st_size
            if expected_size and downloaded_size != expected_size:
                logger.warning(
                    f"Size mismatch: expected {expected_size}, got {downloaded_size}"
                )

            # Rename temp file to final name
            temp_path.rename(output_path)
            logger.info(f"Successfully downloaded: {output_path.name}")
            return True

        except requests.exceptions.RequestException as e:
            logger.warning(f"Attempt {attempt}/{max_retries} failed: {e}")
            if attempt < max_retries:
                wait_time = 2**attempt
                logger.info(f"Retrying in {wait_time} seconds...")
                time.sleep(wait_time)

    logger.error(f"Failed to download after {max_retries} attempts: {url}")
    return False


def download_file_ftp(
    host: str,
    remote_path: str,
    output_path: Path,
    force: bool = False,
    max_retries: int = 3,
    logger: Optional[logging.Logger] = None,
) -> bool:
    """
    Download a file via FTP with progress bar and error handling.

    Args:
        host: FTP host
        remote_path: Path on FTP server
        output_path: Path to save the file
        force: Force re-download even if file exists
        max_retries: Number of retry attempts
        logger: Logger instance

    Returns:
        True if download successful, False otherwise
    """
    if logger is None:
        logger = logging.getLogger("download_references")

    # Check if file already exists
    if output_path.exists() and not force:
        logger.info(f"File already exists, skipping: {output_path.name}")
        return True

    # Create parent directories
    output_path.parent.mkdir(parents=True, exist_ok=True)

    temp_path = output_path.with_suffix(output_path.suffix + ".partial")

    for attempt in range(1, max_retries + 1):
        try:
            with ftplib.FTP(host, timeout=60) as ftp:
                ftp.login()

                # Get file size
                try:
                    file_size = ftp.size(remote_path)
                except Exception:
                    file_size = None

                # Check for resume
                resume_pos = 0
                if temp_path.exists():
                    resume_pos = temp_path.stat().st_size
                    logger.info(f"Resuming from {format_size(resume_pos)}")

                mode = "ab" if resume_pos > 0 else "wb"

                desc = output_path.name[:40]
                with open(temp_path, mode) as f:
                    with tqdm(
                        total=file_size,
                        initial=resume_pos,
                        unit="B",
                        unit_scale=True,
                        unit_divisor=1024,
                        desc=desc,
                        leave=True,
                    ) as pbar:

                        def callback(data: bytes) -> None:
                            f.write(data)
                            pbar.update(len(data))

                        if resume_pos > 0:
                            ftp.retrbinary(
                                f"RETR {remote_path}", callback, rest=resume_pos
                            )
                        else:
                            ftp.retrbinary(f"RETR {remote_path}", callback)

            # Rename temp file to final name
            temp_path.rename(output_path)
            logger.info(f"Successfully downloaded: {output_path.name}")
            return True

        except ftplib.all_errors as e:
            logger.warning(f"Attempt {attempt}/{max_retries} failed: {e}")
            if attempt < max_retries:
                wait_time = 2**attempt
                logger.info(f"Retrying in {wait_time} seconds...")
                time.sleep(wait_time)

    logger.error(f"Failed to download after {max_retries} attempts: {remote_path}")
    return False


def download_1000genomes(
    output_dir: Path, force: bool = False, logger: Optional[logging.Logger] = None
) -> Dict[str, Any]:
    """
    Download 1000 Genomes Project Phase 3 data.

    Args:
        output_dir: Base output directory
        force: Force re-download
        logger: Logger instance

    Returns:
        Dictionary with download statistics
    """
    if logger is None:
        logger = logging.getLogger("download_references")

    info = DATASET_INFO["1000g"]
    dataset_dir = output_dir / info["output_subdir"]
    dataset_dir.mkdir(parents=True, exist_ok=True)

    logger.info("=" * 60)
    logger.info(f"Downloading {info['description']}")
    logger.info(f"Estimated size: ~{info['estimated_size_gb']} GB")
    logger.info(f"Output directory: {dataset_dir}")
    logger.info("=" * 60)

    stats = {"downloaded": 0, "skipped": 0, "failed": 0, "total_size": 0}

    # Parse FTP URL
    parsed = urlparse(info["base_url"])
    host = parsed.netloc
    base_path = parsed.path

    # Download chromosome files
    for filename in info["files"]:
        remote_path = base_path + filename
        output_path = dataset_dir / filename

        success = download_file_ftp(
            host=host,
            remote_path=remote_path,
            output_path=output_path,
            force=force,
            logger=logger,
        )

        if success:
            if output_path.exists():
                stats["total_size"] += output_path.stat().st_size
            stats["downloaded"] += 1
        else:
            stats["failed"] += 1

    # Download panel file
    panel_remote = base_path + info["panel_file"]
    panel_output = dataset_dir / info["panel_file"]
    download_file_ftp(
        host=host,
        remote_path=panel_remote,
        output_path=panel_output,
        force=force,
        logger=logger,
    )

    logger.info(f"\n1000 Genomes download complete:")
    logger.info(f"  Downloaded: {stats['downloaded']} files")
    logger.info(f"  Failed: {stats['failed']} files")
    logger.info(f"  Total size: {format_size(stats['total_size'])}")

    return stats


def download_aadr(
    output_dir: Path, force: bool = False, logger: Optional[logging.Logger] = None
) -> Dict[str, Any]:
    """
    Download Reich Lab Ancient DNA Resource (AADR) data.

    Args:
        output_dir: Base output directory
        force: Force re-download
        logger: Logger instance

    Returns:
        Dictionary with download statistics
    """
    if logger is None:
        logger = logging.getLogger("download_references")

    info = DATASET_INFO["aadr"]
    dataset_dir = output_dir / info["output_subdir"]
    dataset_dir.mkdir(parents=True, exist_ok=True)

    logger.info("=" * 60)
    logger.info(f"Downloading {info['description']}")
    logger.info(f"Estimated size: ~{info['estimated_size_gb']} GB")
    logger.info(f"Output directory: {dataset_dir}")
    logger.info("=" * 60)

    stats = {"downloaded": 0, "skipped": 0, "failed": 0, "total_size": 0}

    for filename in info["files"]:
        url = info["base_url"] + filename
        output_path = dataset_dir / filename

        success = download_file_http(
            url=url, output_path=output_path, force=force, logger=logger
        )

        if success:
            if output_path.exists():
                stats["total_size"] += output_path.stat().st_size
            stats["downloaded"] += 1
        else:
            stats["failed"] += 1

    logger.info(f"\nAADR download complete:")
    logger.info(f"  Downloaded: {stats['downloaded']} files")
    logger.info(f"  Failed: {stats['failed']} files")
    logger.info(f"  Total size: {format_size(stats['total_size'])}")

    return stats


def download_sgdp(
    output_dir: Path, force: bool = False, logger: Optional[logging.Logger] = None
) -> Dict[str, Any]:
    """
    Download Simons Genome Diversity Project (SGDP) data.

    Args:
        output_dir: Base output directory
        force: Force re-download
        logger: Logger instance

    Returns:
        Dictionary with download statistics
    """
    if logger is None:
        logger = logging.getLogger("download_references")

    info = DATASET_INFO["sgdp"]
    dataset_dir = output_dir / info["output_subdir"]
    dataset_dir.mkdir(parents=True, exist_ok=True)

    logger.info("=" * 60)
    logger.info(f"Downloading {info['description']}")
    logger.info(f"Estimated size: ~{info['estimated_size_gb']} GB")
    logger.info(f"Output directory: {dataset_dir}")
    logger.info("=" * 60)

    stats = {"downloaded": 0, "skipped": 0, "failed": 0, "total_size": 0}

    for filename in info["files"]:
        url = info["base_url"] + filename
        output_path = dataset_dir / filename

        success = download_file_http(
            url=url, output_path=output_path, force=force, logger=logger
        )

        if success:
            if output_path.exists():
                stats["total_size"] += output_path.stat().st_size
            stats["downloaded"] += 1
        else:
            stats["failed"] += 1

    logger.info(f"\nSGDP download complete:")
    logger.info(f"  Downloaded: {stats['downloaded']} files")
    logger.info(f"  Failed: {stats['failed']} files")
    logger.info(f"  Total size: {format_size(stats['total_size'])}")

    return stats


def main() -> None:
    """Main entry point for the download script."""
    parser = argparse.ArgumentParser(
        description="Download reference datasets for ancestry analysis.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Datasets and estimated sizes:
  1000g  - 1000 Genomes Project Phase 3 (~15 GB)
  aadr   - Reich Lab Ancient DNA Resource (~2 GB)
  sgdp   - Simons Genome Diversity Project (~500 MB)
  all    - Download all datasets (~17.5 GB total)

Examples:
  %(prog)s --dataset all
  %(prog)s --dataset 1000g
  %(prog)s --dataset aadr --force
  %(prog)s --dataset sgdp --output-dir /path/to/data/
        """,
    )

    parser.add_argument(
        "--dataset",
        type=str,
        required=True,
        choices=["1000g", "aadr", "sgdp", "all"],
        help="Dataset to download: 1000g, aadr, sgdp, or all",
    )

    parser.add_argument(
        "--output-dir",
        type=str,
        default="data/raw",
        help="Output directory for downloaded files (default: data/raw/)",
    )

    parser.add_argument(
        "--force",
        action="store_true",
        help="Force re-download even if files already exist",
    )

    parser.add_argument(
        "--log-dir",
        type=str,
        default="logs",
        help="Directory for log files (default: logs/)",
    )

    args = parser.parse_args()

    # Set up paths
    output_dir = Path(args.output_dir).resolve()
    log_dir = Path(args.log_dir).resolve()

    # Set up logging
    logger = setup_logging(log_dir)

    # Print banner
    logger.info("=" * 60)
    logger.info("Ancestry Analysis Pipeline - Reference Data Downloader")
    logger.info("=" * 60)
    logger.info(f"Dataset: {args.dataset}")
    logger.info(f"Output directory: {output_dir}")
    logger.info(f"Force re-download: {args.force}")
    logger.info("")

    start_time = time.time()
    total_stats = {"downloaded": 0, "failed": 0, "total_size": 0}

    # Download requested datasets
    datasets_to_download: List[str] = []
    if args.dataset == "all":
        datasets_to_download = ["1000g", "aadr", "sgdp"]
    else:
        datasets_to_download = [args.dataset]

    for dataset in datasets_to_download:
        try:
            if dataset == "1000g":
                stats = download_1000genomes(output_dir, args.force, logger)
            elif dataset == "aadr":
                stats = download_aadr(output_dir, args.force, logger)
            elif dataset == "sgdp":
                stats = download_sgdp(output_dir, args.force, logger)
            else:
                logger.error(f"Unknown dataset: {dataset}")
                continue

            total_stats["downloaded"] += stats["downloaded"]
            total_stats["failed"] += stats["failed"]
            total_stats["total_size"] += stats["total_size"]

        except KeyboardInterrupt:
            logger.warning("\nDownload interrupted by user")
            sys.exit(1)
        except Exception as e:
            logger.error(f"Error downloading {dataset}: {e}")
            total_stats["failed"] += 1

    # Print summary
    elapsed_time = time.time() - start_time
    hours, remainder = divmod(int(elapsed_time), 3600)
    minutes, seconds = divmod(remainder, 60)

    logger.info("")
    logger.info("=" * 60)
    logger.info("Download Summary")
    logger.info("=" * 60)
    logger.info(f"Total files downloaded: {total_stats['downloaded']}")
    logger.info(f"Total files failed: {total_stats['failed']}")
    logger.info(f"Total size: {format_size(total_stats['total_size'])}")
    logger.info(f"Time elapsed: {hours:02d}:{minutes:02d}:{seconds:02d}")
    logger.info(f"Log file: {log_dir / 'download.log'}")
    logger.info("=" * 60)

    if total_stats["failed"] > 0:
        logger.warning(
            f"Some downloads failed. Check the log file for details."
        )
        sys.exit(1)


if __name__ == "__main__":
    main()
