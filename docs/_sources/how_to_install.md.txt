# HTMS_Amber: High-Throughput Mutational Scanning Pipeline with mmpbsa.py 
---
### Summary
This package aims to improve upon the base Alanine scanning method within [Amber](https://ambermd.org/) by offering a rapid and easily accessible computational tool for calculating changes in Binding Free Energy ($\Delta\Delta G_{bind}$) via mmpbsa.py for in silico mutations. 

The package supports all amino acid mutations, including single-point and multi-point mutations. Furthermore, a case study involving SARS-CoV-2 variants and the human ACE2 receptor is provided in the examples section.

For more details on running simulations, see the Usage Section below.
---

# Installation

### 1. Prerequisites
You will need:
- **Conda** or **Mamba**: If you do not have it, see [Miniforge](https://github.com/conda-forge/miniforge).
- [**AMBER**](https://ambermd.org/): Required for Alanine Scanning and MMPBSA calculations.
- [**MODELLER**](https://salilab.org/modeller/): Required for multiple point mutations and specific single point mutations (see {ref}`non-alanine-mutations`)
<!-- (see [Mutation Details](#-non-alanine-mutations)).  -->
### 2. Clone the Repository
Download the project to your local machine using one of the following methods:

**Via GitHub Desktop:**
1. Click the green **Code** button at the top of this page.
2. Select **Open with GitHub Desktop**.

**Via Command Line:**
```bash
git clone #TODO: UPDATE
cd HTMS_Amber
```

### 3. Create and activate the conda env
```bash 
conda env create -f htms_conda_env.yml
conda activate amber_muts
```

### Notes on **MODELLER** installation
To make full use of the pipline a reistered MODELLER install is required see [modeller registration](https://salilab.org/modeller/registration.html )

#### Resource Optimization
1. Use GPU resources for the production MD runs.
2. Scale back to minimal CPU and memory for the MMPBSA.py analysis portion.




