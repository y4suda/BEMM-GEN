# Celluar Environment Mimicking Model Geneletaor (CEMM-GEN)

## Requirements 


## Files  


## How to use
Now that we have prepared files for the M2 receptor, agonist iperoxo and initial water placement, we need to build the membrane simulation box.

cd ../membrane_build

We run PACKMOL-Memgen using the "m2_prep.pdb" file as input:

packmol-memgen --pdb ../system_pdb/m2_prep.pdb --lipids POPC:CHL1 --ratio 9:1 --preoriented --salt --salt_c Na+ --saltcon 0.15 --dist 10 --dist_wat 15 --notprotonate --nottrim

Take care to understand each of the flags. Here we ask for a mixed POPC/CHOL membrane, but you may want plain POPC, or a different lipid composition. You can run packmol-memgen --help to find out about each flag. Briefly:

"--pdb ../system_pdb/m2_prep.pdb": input receptor PDB
"--lipids POPC:CHL1": build a mixed POPC/CHOL membrane
"--ratio 9:1": POPC/CHOL ratio of 9:1
"--preoriented": the OPM coordinates are already orientated along the z-axis, for membrane building
"--salt --salt_c Na+ --saltcon 0.15": Add 0.15 M of NaCl salt to the water layer
"--dist 10": minimum maxmin value for x, y, z to box boundary
"--dist_wat 15": water layer thickness of 15 A
"--notprotonate --nottrim": do not process input receptor PDB file (since we have already prepared this)



## Author 
Takunori Yasuda, Rikuri Morita, Yasuteru Shigeta, Ryuhei Harada.  
Center for Computational Sciences, University of Tsukuba  
takunoriyasuda@gmail.com
