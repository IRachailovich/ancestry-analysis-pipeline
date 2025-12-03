# Data Sources and Download Instructions

This document provides information about where to download reference datasets for ancestry analysis and how to properly cite them.

## Table of Contents
1. [Reich Lab Ancient DNA Datasets](#1-reich-lab-ancient-dna-datasets)
2. [Modern Reference Populations](#2-modern-reference-populations)
3. [SNP Panels and Coordinates](#3-snp-panels-and-coordinates)
4. [Download Automation](#4-download-automation)
5. [Data Organization](#5-data-organization)
6. [Citation Guidelines](#6-citation-guidelines)

---

## 1. Reich Lab Ancient DNA Datasets

### Human Origins Dataset
The Human Origins array is the standard SNP panel used in ancient DNA studies.

**Download**: https://reich.hms.harvard.edu/allen-ancient-dna-resource-aadr-downloadable-genotypes-present-day-and-ancient-dna-data

**Contents**:
- ~600,000 SNPs
- Modern reference populations from diverse global origins
- Ancient DNA samples spanning 45,000 years
- Pre-computed genetic data in EIGENSTRAT and PLINK formats

**Files to download**:
- v54.1.p1_1240K_public.tar (Full 1240k dataset)
- v54.1.p1_HO_public.tar (Human Origins subset ~600k SNPs)

**Citation**:
- Lazaridis et al. (2014). "Ancient human genomes suggest three ancestral populations for present-day Europeans." Nature 513, 409–413.
- Patterson et al. (2012). "Ancient admixture in human history." Genetics 192, 1065–1093.

### Allen Ancient DNA Resource (AADR)
Comprehensive database of published ancient genomes.

**Download**: Same link as above (updated regularly)

**Citation**: Cite the AADR and the original publications for specific samples you use.

---

## 2. Modern Reference Populations

### 1000 Genomes Project (Phase 3)
High-coverage whole genome sequences from 26 populations worldwide.

**Download**: 
- FTP: ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/
- Alternative: https://www.internationalgenome.org/data

**Populations**:
- African (YRI, LWK, GWD, MSL, ESN, ASW, ACB)
- European (CEU, TSI, FIN, GBR, IBS)
- East Asian (CHB, JPT, CHS, CDX, KHV)
- South Asian (GIH, PJL, BEB, STU, ITU)
- American (MXL, PUR, CLM, PEL)

**Format**: VCF files by chromosome

**Citation**:
- 1000 Genomes Project Consortium (2015). "A global reference for human genetic variation." Nature 526, 68–74.

### Simons Genome Diversity Project (SGDP)
High-coverage genomes from 300 individuals representing 142 diverse populations.

**Download**: https://sharehost.hms.harvard.edu/genetics/reich_lab/sgdp/

**Citation**:
- Mallick et al. (2016). "The Simons Genome Diversity Project: 300 genomes from 142 diverse populations." Nature 538, 201–206.

---

## Next Steps

After downloading data:
1. Proceed to 03_workflow.md for analysis instructions
2. Convert your 23andMe data to PLINK format
3. Merge with reference populations
4. Run quality control
