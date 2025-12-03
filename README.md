# Ancestry Analysis Pipeline

A comprehensive population genetics pipeline for analyzing modern SNP array data (23andMe, AncestryDNA) and ancient DNA to infer ancestry, admixture proportions, and population history.

## 🎯 Project Goals

This pipeline enables you to:

1. **Infer global ancestry proportions** from modern SNP array data using PCA and ADMIXTURE
2. **Estimate ancient ancestry components** (Yamnaya, Early European Farmers, Western Hunter-Gatherers) in modern and ancient samples using qpAdm
3. **Perform local ancestry inference** along chromosomes to identify ancestry-specific genomic segments using RFMix
4. **Compare populations and build admixture graphs** to model population history using TreeMix and qpGraph

## 📊 Data Sources

This project integrates:
- **Your personal genetic data**: 23andMe raw data
- **Ancient DNA**: Reich Lab datasets (Human Origins, 1240k SNP panel)
- **Modern reference populations**: 1000 Genomes Project, Simons Genome Diversity Project (SGDP)

## 🛠️ Tech Stack

- **Python 3.10+**: Main pipeline language
- **Population genetics tools**: PLINK, ADMIXTURE, EIGENSOFT (smartpca), ADMIXTOOLS (qpAdm/qpDstat/qpGraph), RFMix, TreeMix
- **Analysis & visualization**: pandas, numpy, matplotlib, seaborn, scikit-allel
- **Environment**: Windows with WSL2 (Windows Subsystem for Linux)

## 📁 Repository Structure

```
ancestry-analysis-pipeline/
├── README.md                          # This file
├── docs/                              # Documentation
│   ├── 01_setup.md                   # Software installation (WSL2, tools)
│   ├── 02_data_sources.md            # Public dataset download instructions
│   ├── 03_workflow.md                # Step-by-step analysis workflow
│   └── 04_interpretation.md          # Results interpretation guide
├── data/                              # Data directory (gitignored)
│   ├── raw/                          # Raw input data
│   │   ├── 23andme/                  # Your 23andMe data
│   │   ├── reich_lab/                # Ancient DNA datasets
│   │   └── modern_references/        # 1000 Genomes, SGDP
│   ├── processed/                    # QC'd, filtered, merged datasets
│   └── reference_panels/             # Curated reference populations
├── scripts/                           # Analysis scripts
│   ├── 00_setup/                     # Setup and data download
│   ├── 01_preprocessing/             # Data QC and format conversion
│   ├── 02_global_ancestry/           # PCA and ADMIXTURE
│   ├── 03_ancient_ancestry/          # qpAdm analysis
│   ├── 04_local_ancestry/            # RFMix local ancestry
│   ├── 05_admixture_graphs/          # TreeMix and qpGraph
│   └── utils/                        # Shared utility functions
├── results/                           # Analysis outputs (gitignored)
│   ├── pca/
│   ├── admixture/
│   ├── qpadm/
│   ├── local_ancestry/
│   └── admixture_graphs/
├── notebooks/                         # Jupyter notebooks for exploration
├── environment.yml                    # Conda environment specification
├── requirements.txt                   # Python package dependencies
└── .gitignore                        # Git ignore rules (excludes data files)
```

## 🚀 Quick Start

### Prerequisites (Windows)

1. **Install WSL2** (Windows Subsystem for Linux):
   ```powershell
   wsl --install
   ```
   Restart your computer after installation.

2. **Install Miniconda in WSL2**:
   ```bash
   wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
   bash Miniconda3-latest-Linux-x86_64.sh
   ```

### Installation

1. **Clone this repository**:
   ```bash
   git clone https://github.com/IRachailovich/ancestry-analysis-pipeline.git
   cd ancestry-analysis-pipeline
   ```

2. **Create the conda environment**:
   ```bash
   conda env create -f environment.yml
   conda activate ancestry-analysis
   ```

3. **Install additional tools**:
   ```bash
   bash scripts/00_setup/install_tools.sh
   ```

4. **Download reference datasets**:
   ```bash
   python scripts/00_setup/download_references.py
   ```

### Running the Analysis

1. **Place your 23andMe data** in `data/raw/23andme/`

2. **Convert and QC your data**:
   ```bash
   python scripts/01_preprocessing/convert_23andme.py
   python scripts/01_preprocessing/qc_filtering.py
   ```

3. **Run global ancestry analysis**:
   ```bash
   python scripts/02_global_ancestry/run_pca.py
   python scripts/02_global_ancestry/run_admixture.py
   ```

4. **Estimate ancient ancestry proportions**:
   ```bash
   python scripts/03_ancient_ancestry/prepare_qpadm.py
   bash scripts/03_ancient_ancestry/run_qpadm.sh
   ```

5. **Explore results** in Jupyter notebooks:
   ```bash
   jupyter notebook notebooks/
   ```

## 📚 Documentation

- **[Setup Guide](docs/01_setup.md)**: Detailed installation instructions for WSL2 and all tools
- **[Data Sources](docs/02_data_sources.md)**: Where to download reference datasets and how to cite them
- **[Workflow Guide](docs/03_workflow.md)**: Step-by-step walkthrough of each analysis
- **[Interpretation Guide](docs/04_interpretation.md)**: How to read and interpret your results

## 🔬 Methods & Tools

### Global Ancestry
- **PCA**: EIGENSOFT (smartpca) for population structure visualization
- **ADMIXTURE**: Model-based clustering for ancestry proportions (K=2 to K=15)

### Ancient Ancestry
- **qpAdm** (ADMIXTOOLS): Estimate proportions from ancient source populations
- **D-statistics / f-statistics**: Test admixture hypotheses

### Local Ancestry
- **RFMix**: Chromosome painting and local ancestry inference

### Admixture Graphs
- **TreeMix**: Maximum likelihood phylogenetic trees with migration edges
- **qpGraph**: Admixture graph fitting and testing

## 📖 Key References

- **ADMIXTURE**: Alexander, Novembre & Lange (2009). *Genome Research*
- **qpAdm**: Patterson et al. (2012). *Genetics*; Haak et al. (2015). *Nature*
- **RFMix**: Maples et al. (2013). *AJHG*
- **Human Origins dataset**: Lazaridis et al. (2014). *Nature*; Patterson et al. (2012)

See [docs/04_interpretation.md](docs/04_interpretation.md) for full bibliography.

## ⚠️ Privacy & Data Security

- **All genetic data files are gitignored** and will never be committed to GitHub
- Your 23andMe data stays local on your machine
- Only scripts and documentation are version controlled
- When sharing results, always anonymize sample IDs

## 📝 License

This project is for personal research and educational purposes. Public datasets have their own licenses (see [docs/02_data_sources.md](docs/02_data_sources.md)).

## 🤝 Contributing

This is a personal research project, but suggestions and improvements are welcome via issues.

## 📧 Contact

GitHub: [@IRachailovich](https://github.com/IRachailovich)

---

**Ready to explore your ancestry? Start with [docs/01_setup.md](docs/01_setup.md)!**