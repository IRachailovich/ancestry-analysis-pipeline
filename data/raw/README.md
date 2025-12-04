# Reference Datasets for Ancestry Analysis

This directory contains reference datasets used for ancestry analysis in this pipeline.

## Dataset Overview

### 1. 1000 Genomes Project Phase 3 (`1000genomes/`)

High-coverage whole genome sequences from 26 populations worldwide.

**Source**: ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/

**Contents**:
- VCF files for chromosomes 1-22
- Population panel file (`integrated_call_samples_v3.20130502.ALL.panel`)

**Populations included**:
- **African**: YRI, LWK, GWD, MSL, ESN, ASW, ACB
- **European**: CEU, TSI, FIN, GBR, IBS
- **East Asian**: CHB, JPT, CHS, CDX, KHV
- **South Asian**: GIH, PJL, BEB, STU, ITU
- **American**: MXL, PUR, CLM, PEL

**Citation**:
> 1000 Genomes Project Consortium (2015). "A global reference for human genetic variation." *Nature* 526, 68–74. https://doi.org/10.1038/nature15393

---

### 2. Allen Ancient DNA Resource (AADR) (`aadr/`)

Comprehensive database of published ancient genomes from the Reich Lab.

**Source**: https://reichdata.hms.harvard.edu/pub/datasets/amh_repo/curated_releases/

**Contents**:
- `v54.1_1240K_public.geno` - Genotype data
- `v54.1_1240K_public.snp` - SNP information
- `v54.1_1240K_public.ind` - Individual information

**Key ancient populations**:
- Yamnaya (Steppe pastoralists)
- Anatolia_N (Early Anatolian farmers)
- WHG (Western Hunter-Gatherers)
- EHG (Eastern Hunter-Gatherers)
- CHG (Caucasus Hunter-Gatherers)

**Format**: EIGENSTRAT

**Citations**:
> Lazaridis et al. (2014). "Ancient human genomes suggest three ancestral populations for present-day Europeans." *Nature* 513, 409–413.

> Patterson et al. (2012). "Ancient admixture in human history." *Genetics* 192, 1065–1093.

---

### 3. Simons Genome Diversity Project (SGDP) (`sgdp/`)

High-coverage genomes from 300 individuals representing 142 diverse populations worldwide.

**Source**: https://reichdata.hms.harvard.edu/pub/datasets/sgdp/

**Contents**:
- `cteam_extended.v4.geno` - Genotype data
- `cteam_extended.v4.snp` - SNP information
- `cteam_extended.v4.ind` - Individual information

**Format**: EIGENSTRAT

**Citation**:
> Mallick et al. (2016). "The Simons Genome Diversity Project: 300 genomes from 142 diverse populations." *Nature* 538, 201–206. https://doi.org/10.1038/nature18964

---

## Download Instructions

Use the download script to automatically fetch these datasets:

```bash
# Download all datasets
python scripts/00_setup/download_references.py --dataset all

# Download individual datasets
python scripts/00_setup/download_references.py --dataset 1000g
python scripts/00_setup/download_references.py --dataset aadr
python scripts/00_setup/download_references.py --dataset sgdp

# Force re-download
python scripts/00_setup/download_references.py --dataset all --force
```

## Data Formats

### EIGENSTRAT Format
Used by AADR and SGDP datasets:
- `.geno` - Genotype matrix (0/1/2/9 encoding)
- `.snp` - SNP information (ID, chromosome, position, alleles)
- `.ind` - Individual information (ID, sex, population)

### VCF Format
Used by 1000 Genomes:
- Standard Variant Call Format
- Compressed with gzip (`.vcf.gz`)

### PLINK Format
May need to convert VCF to PLINK binary format:
- `.bed` - Binary genotype data
- `.bim` - SNP information
- `.fam` - Individual information

## Storage Requirements

| Dataset | Estimated Size |
|---------|---------------|
| 1000 Genomes | ~15 GB |
| AADR | ~2 GB |
| SGDP | ~500 MB |
| **Total** | **~17.5 GB** |

## Notes

- Large files are not committed to git (see `.gitignore`)
- Downloads support resuming interrupted transfers
- Checksums are verified when available
- See `logs/download.log` for download history

## Privacy Notice

Reference datasets are public and do not contain personally identifiable information.
Your personal genetic data (23andMe, etc.) should be stored in `data/raw/23andme/`
and is excluded from version control.
