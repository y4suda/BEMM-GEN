<p align="center"><br><img src="./cemm-gen_logo_horizontal.png" width="512px"><br><br>
Cellular Environment Mimicking Model GENerator <br><br><br></p>

# Welcome to CEMM-GEN
<p align="center"><br><img src="./cemm-gen_main.png" width="800px"></p>


## What is CEMM-GEN?
<img src="./cemm-gen_scheme.png" align="right" width="400px">
CEMM-GEN is a new tool aiming to generate cellular environment mimicking model for molecular dynamics simulation. CEMM-GEN generates spherical and cylindrical models with user-specified chemical properties, allowing the integration of arbitrary protein structures into the generated models. Consequently, the model and protein complex structures, along with the corresponding parameter files for MD simulations (AMBER or GROMACS), are provided as output files. 
<br><br><br><br><br><br><br><br><br>

## Requirements
Python >= 3.8, conda, openbabel, psi4, resp, ambertools, scikit-learn, rdkit, parmed 

## Installation 
### For Linux/MacOS
- Follow the instructions below to install CEMM-GEN.
### For Windows
- Install Windows Subsystem for Linux (WSL), then follow the instructions below to install CEMM-GEN on WSL.


Using a Conda environment (such as Miniforge, or Miniconda) is recommended for CEMM-GEN. If you use only the pip package manager, please adjust the commands accordingly.
```sh
conda create -n cemm-gen-env
conda activate cemm-gen-env
conda install -c conda-forge openbabel psi4 resp ambertools
pip install cemm-gen
```
## Basic Usage
```sh:available_sub-commands
# Make a cylindrical model
cemm-gen cylinder --proteinseq GASGASGASGAS --proteinSS HHHHHHHHHHHH --resnames MTY:HYD --composition 1:2.5

# Make a spherical model
cemm-gen sphere --proteinpdb protein.pdb --resnames MTY:HYD --composition 0.3:0.7

# Make parameters for a new residue
cemm-gen makeparam --smiles CCC --resname MTY --description "Methyl group"
```

**Please see the documentation [English](https://github.com/y4suda/CEMM-GEN/blob/main/tutorial_en.md) / [日本語](https://github.com/y4suda/CEMM-GEN/blob/main/tutorial_ja.md)**

## Cite Us
T. Yasuda, R. Morita, Y. Shigeta and R. Harada. "Cellular Environment Mimicking Model GENerator: A tool for generating a cellular environment mimicking model."*XXXX*,2024,XXXX,[doi](https://XXX)

## Authors
Takunori Yasuda, Doctoral Program in Biology, University of Tsukuba
Rikuri Morita, Center for Computational Sciences, University of Tsukuba  
Yasuteru Shigeta, Center for Computational Sciences, University of Tsukuba  
Ryuhei Harada, Center for Computational Sciences, University of Tsukuba  
yasuda.takunori.tkb_gb@u.tsukuba.ac.jp 
