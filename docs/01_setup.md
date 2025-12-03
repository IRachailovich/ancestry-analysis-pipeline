# Setup Guide for Windows WSL2 Installation

## Introduction
This document provides a comprehensive setup guide for running the Ancestry Analysis Pipeline on Windows using Windows Subsystem for Linux (WSL2).

## WSL2 Installation
1. Enable Windows Subsystem for Linux from the Windows features.
2. Install WSL2 by following the official Microsoft guide.

## Miniconda Installation
1. Download the Miniconda installer for Linux.
2. Run the installer and follow the prompts to complete the installation.

## Conda Environment Setup
1. Create a new conda environment using the command:
   ```bash
   conda create --name myenv python=3.8
   ```
2. Activate the conda environment:
   ```bash
   conda activate myenv
   ```

## Installation of EIGENSOFT (smartpca)
1. Follow the EIGENSOFT installation guide to set up smartpca within your conda environment.

## Installation of ADMIXTOOLS (qpAdm/qpDstat/qpGraph)
1. Clone the ADMIXTOOLS repository:
   ```bash
   git clone <repository-url>
   ```
2. Follow the instructions to compile and install tools:
   ```bash
   make
   ```

## RFMix Installation
1. Follow the official RFMix installation instructions according to your environment.

## TreeMix Installation
1. Follow the installation guide for TreeMix on the official GitHub page.

## Verification Commands
- After installations, verify each tool with the following commands:
   ```bash
   smartpca -h
   qpAdm -h
   RFMix -h
   TreeMix -h
   ```

## Troubleshooting Section
### WSL2 Issues
- Common issues and fixes for WSL2 setup:
  - Ensure virtualization is enabled in BIOS.
  - Check for Windows updates.

### Conda Issues
- Solutions for common conda problems:
  - Use `conda clean --all` to clear cache.

### Compilation Errors
- Steps to resolve compilation errors when building ADMIXTOOLS.

### PATH Problems
- Ensure that the path to binaries is correctly set in your `.bashrc`.

### Memory Allocation
- For WSL2 memory allocation, modify your `.wslconfig` file:
   ```
   [wsl2]
   memory=4GB  # Set the amount of memory
   processors=2 # Set number of processors
   ```

---
This guide is intended to provide a comprehensive approach to set up the ancestry analysis tools effectively in a WSL2 environment.