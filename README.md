# Ribotunnel 

*Jupyter notebooks and Python code used in the paper: "Evolution of the ribosomal exit tunnel through the eyes of the nascent chain"*

## Overview
This repository contains the scripts and notebooks required to reproduce the analysis and build the models described in the paper. 

## Repository Structure

### Python Modules
The core Python functions and scripts are located in the `ribosome_analysis` directory.

### Jupyter Notebooks
- **`Methodology.ipynb`** -> Methodology used to build the ribosome tunnel model for simulations and set up the MD simulations.
- **`PTC.ipynb`** -> Used to calculate the Peptidyl Transferase Center (PTC) RMSD across many different ribosome structures.
- **`Trajectory_analysis.ipynb`** -> Step-by-step pipeline for the MD trajectory analysis.
- **`tunnel_analysis_cytoplasm.ipynb`** -> Tools for generating tunnels from MD simulations and carrying out the complete structural analysis.

## Requirements
VMD, ChimeraX, GROMACS

## Reference
If you use this code or build upon this work, please cite:

> Wlodarski T., "Evolution of the ribosomal exit tunnel through the eyes of the nascent chain", bioRxiv 2025.  
> [Read the preprint on bioRxiv](https://www.biorxiv.org/content/10.64898/2025.12.25.696519v2)
