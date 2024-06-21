<p align="center"><br><img src="./cemm-gen_logo_horizontal.png" width="512px"><br><br>
Cellular Environment Mimicking Model GENerator  ver. 2024.4.1<br><br><br></p>

# Welcome to CEMM-GEN

## What is CEMM-GEN?
CEMM-GEN is a new tool aiming to generate cellular environment mimicking model for molecular dynamics simulation. CEMM-GEN provides a configuration with any modification group added to the spherical or cylinder-like interior. 
##
___
## Instllation 
```sh
conda create -n cemm-gen-env
conda activate cemm-gen-env
conda install -c conda-forge openbabel psi4 resp ambertools
pip install -e {CEMM-GEN_dir}
```
## Basic Usage
```sh:available_sub-commands
# Make a cylindrical model
cemm-gen cylinder

# Make a spherical model
cemm-gen sphere

# Make parameters for a new residue
cemm-gen makeparam
```

**Please see the documentation [here](https://github.com/y4suda/CEMM-GEN/blob/main/tutorial_en.md)**
___

## Cite Us
T. Yasuda, R. Morita, Y. Shigeta and R. Harada. "Cellular Environment Mimicking Model GENerator: A tool for generating a cellular environment mimicking model."*XXXX*,2024,XXXX,[doi](https://XXX)

## Authors
Takunori Yasuda, Doctoral Program in Biology  
Rikuri Morita, Center for Computational Sciences, University of Tsukuba  
Yasuteru Shigeta, Center for Computational Sciences, University of Tsukuba  
Ryuhei Harada, Center for Computational Sciences, University of Tsukuba  
yasuda.takunori.tkb_gb@u.tsukuba.ac.jp 
