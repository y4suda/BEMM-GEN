import argparse
import subprocess
import os
import shutil

import numpy as np
from scipy.spatial import distance
import parmed

from . import utils
from . import make_param


def make_MD_input(args: argparse.Namespace):

    # remove hydrogens in the protein.pdb
    shutil.copyfile("protein.pdb", "protein_H.pdb")
    with open("protein.pdb", "w") as f:
        subprocess.run("reduce -Trim protein_H.pdb", shell=True, check=False, stdout=f, stderr=subprocess.DEVNULL)
    os.remove("protein_H.pdb")

    _covert_amber(args)
    _restraint_amber(args)

    _convert_GROMACS(args)
    _restraint_GROMACS(args)
    _remove_overlap(f"{args.output_prefix}.gro")

    return True

def _covert_amber(args: argparse.Namespace):
    leap_command = ""
    leap_command += f"source leaprc.protein.{args.forcefield_protein}\n"
    leap_command += f"source leaprc.water.tip3p\n"
    leap_command += f"source leaprc.{args.forcefield_model}\n"
    
    available_fr_list=make_param._get_params()    
    for resname in args.resnames.split(":"):
        if resname not in available_fr_list:
            raise ValueError(f"{resname} is not available")
        if os.path.getsize(f"{available_fr_list[resname]}/{resname}.lib") != 0:
            leap_command += f"loadamberparams {available_fr_list[resname]}/{resname}.frcmod\n"
            leap_command += f"loadoff {available_fr_list[resname]}/{resname}.lib\n"

    if args.command == "cylinder":
        leap_command += f"model = loadpdb cylinder.pdb\n"
    elif args.command == "sphere":
        leap_command += f"model = loadpdb sphere.pdb\n"

    if args.proteinpdb is not None:
        leap_command += f"protein = loadpdb {args.proteinpdb}\n"
        # leap_command += f"system = protein\n"
        leap_command += f"system = combine {{ protein model }}\n"
    else:
        leap_command += f"system = model\n"

    leap_command += f"solvateBox system TIP3PBOX 10 \n"
    leap_command += f"addions system Na+ 0\n"
    leap_command += f"addions system Cl- 0\n"
    leap_command += f"charge system\n"

    leap_command += f"saveamberparm system {args.output_prefix}.prmtop {args.output_prefix}.inpcrd\n"
    leap_command += f"savepdb system {args.output_prefix}.pdb\n"
    leap_command += f"quit\n"

    with open("leap.in", "w") as f:
        f.write(leap_command)
    
    exitcode = subprocess.run("tleap -f leap.in", shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    if exitcode.returncode != 0:
        utils.print_error("Failed to convert pdb file to amber format in tleap. See leap.log for details.")
        exit(1)
    else:
        utils.print_info(f"Amber files are created and saved. ({args.output_prefix}.prmtop and {args.output_prefix}.inpcrd)")

    return True

def _convert_GROMACS(args: argparse.Namespace):

    parm = parmed.load_file(f"{args.output_prefix}.prmtop", f"{args.output_prefix}.inpcrd")
    parm.save(f"{args.output_prefix}.top", format="gromacs", overwrite=True)
    parm.save(f"{args.output_prefix}.gro", overwrite=True)

    utils.print_info(f"GROMACS files are created and saved. ({args.output_prefix}.top and {args.output_prefix}.gro)")
    return True

def _restraint_amber(args: argparse.Namespace):
    utils.print_info("Use the following options in .in files for the Amber simulation.")
    model_start = 0
    model_end = 0
    protein_end = 0
    water_start = 0
    is_protein = True
    with open(f"{args.output_prefix}.pdb", "r") as f:
        for line in f:
            if line.startswith("ATOM"):
                if is_protein:
                    protein_end = int(line[22:26])
                if water_start == 0 and line[17:20] == "WAT":
                    water_start = int(line[22:26])
            if line.startswith("TER"):
                is_protein = False

    restraint_force = [10000, 10000, 10000]
    print("\n------ amber_restraint.in ------")
    print(f"ntr=1,\nrestraintmask=':{protein_end+1}-{water_start-1} & @C,C1',\nrestraint_wt={restraint_force[0]}")
    print("---------------------------------\n")
    with open("amber_restraint.in", "w") as writer:
        writer.write(f"ntr=1,\nrestraintmask=':{protein_end+1}-{water_start-1} & @C,C1',\nrestraint_wt={restraint_force[0]}")

    return True

def _restraint_GROMACS(args: argparse.Namespace):

    restraint_force = [10000, 10000, 10000]
    protein_restraint_force = [1000, 1000, 1000]
    shutil.move(f"{args.output_prefix}.top", f"{args.output_prefix}_norestraint.top")
    
    residue_list = args.resnames.split(":")
    with open(f"{args.output_prefix}_norestraint.top", "r") as f:
        with open(f"{args.output_prefix}.top", "w") as g:
            is_restraint_molecule = False
            is_restraint_protein = None
            restraint_protein_index_list = []
            for line in f:
                if line.startswith("[ moleculetype ]"):

                    if is_restraint_protein == None:
                        is_restraint_protein = True
                    elif is_restraint_protein == True:
                        is_restraint_protein = False
                        g.write("#ifdef PROTEINPOSRES\n")
                        g.write("[ position_restraints ]\n")
                        for i in restraint_protein_index_list:
                            g.write(f"{i: >5}    1    {protein_restraint_force[0]} {protein_restraint_force[1]} {protein_restraint_force[2]}\n")
                        g.write("#endif\n")
                        g.write("\n")
                    else:
                        pass


                    # Add position restraints of previous molecule
                    if is_restraint_molecule == True:
                        if "system" in molecule_name:
                            # restraint for ACE
                            g.write("[ position_restraints ]\n")
                            g.write(f"{2: >5}    1    {restraint_force[0]} {restraint_force[1]} {restraint_force[2]}\n")
                            g.write(f"{5: >5}    1    {restraint_force[0]} {restraint_force[1]} {restraint_force[2]}\n")
                            g.write("\n")
                        else:
                            g.write("[ position_restraints ]\n")
                            g.write(f"{1: >5}    1    {restraint_force[0]} {restraint_force[1]} {restraint_force[2]}\n")
                            g.write(f"{2: >5}    1    {restraint_force[0]} {restraint_force[1]} {restraint_force[2]}\n")
                            g.write("\n")
                        is_restraint_molecule = False


                    g.write(line)
                    g.write(f.readline())
                    molecule_name_line = f.readline()
                    g.write(molecule_name_line)
                    molecule_name = molecule_name_line.split()[0]

                    # ligand-like residues
                    if molecule_name in residue_list:
                        is_restraint_molecule = True

                    # peptide-like residues
                    if molecule_name != "system1" and "system" in molecule_name:
                        is_restraint_molecule = True
                    
                    continue
                    
                if is_restraint_protein == True:
                    split_line = line.split()
                    if len(split_line) == 11 and split_line[4] == "CA":
                        restraint_protein_index_list.append(int(split_line[0]))

                g.write(line)
    utils.print_info("Position restraints were added in topology file for the GROMACS simulation.")

    return True

def _remove_overlap(filename: str):
    gro_file_name = filename
    shutil.copy(gro_file_name, gro_file_name.replace(".gro","_orig.gro"))

    #座標リストの作成
    atom_coordinate = []
    with open(gro_file_name,"r") as f:
        for x,line in enumerate(f):
            if x > 1 and len(line.split()) != 3:
                atom_coordinate.append(line.split()[-3:])
    atom_coordinate_float = []
    for i in atom_coordinate:
        atom_coordinate_float.append([float(j) for j in i])


    block_size = min(len(atom_coordinate_float), 10000)
    is_overlap =[[True] *((len(atom_coordinate_float)-1)//block_size+1) for i in range((len(atom_coordinate_float)-1)//block_size+1)]

    utils.print_info("Removing atomic clashes...")

    removing_cycle = 10
    for cycle in range(removing_cycle):
        print(f"Starting cycle {cycle} ... ", end="")
        current_overlap_num = 0
        #距離行列の作成
        for block_A in range(0, (len(atom_coordinate_float)-1)//block_size+1):
            for block_B in range(0, (len(atom_coordinate_float)-1)//block_size+1):
                dist = distance.cdist(atom_coordinate_float[block_A*block_size:(block_A+1)*block_size], atom_coordinate_float[block_B*block_size:(block_B+1)*block_size], 'euclidean')
                #オーバーラップの検出
                cutoff = 0.05
                overlap = np.where((dist > 0)&(dist < cutoff))
                overlap_list = list(set(zip(*overlap)))
                current_overlap_num += len(overlap_list)
                overlap_removed = set([str(sorted(i)[0]+block_A*block_size) + "|" + str(sorted(i)[1]+block_B*block_size) for i in overlap_list])
                overlap_removed_list = [[ int(j) for j in (i.split("|"))] for i in overlap_removed]

                if len(overlap_removed_list) != 0:
                    is_overlap[block_A][block_B] = True
                else:
                    is_overlap[block_A][block_B] = False

                for j in overlap_removed_list:
                    atom_a = np.array(atom_coordinate_float[j[0]])
                    atom_b = np.array(atom_coordinate_float[j[1]])
                    vector = atom_a - atom_b
                    distance_a_b = np.linalg.norm(vector)
                    if distance_a_b != 0:
                        atom_b_n = atom_b - vector/distance_a_b * (cutoff/2)
                    else:
                        atom_b_n = atom_b + np.random.rand(3) * (cutoff/2)
                    atom_coordinate_float[j[1]] = atom_b_n
        print(f" still {current_overlap_num} overlaps.")
        if not any(np.array(is_overlap).flatten()):
            utils.print_info("Complete removing atomic crashes. There is no overlap.") 
            break
        if cycle == removing_cycle-1:
            utils.print_warning(f"There are still overlaps after {removing_cycle} cycles. Please check the output file carefully.")
            utils.print_info("Changing options \"--min-distance\", and \"--padding-radius\" to larger values might help to avoid overlaps.")

    with open(gro_file_name, "w") as writer:
        with open(gro_file_name.replace(".gro","_orig.gro"), "r") as f:
            for x, line in enumerate(f):
                if x > 1 and len(line.split()) != 3 and "WAT" not in line and "SOL" not in line:
                    line = line[:20] +f'{(atom_coordinate_float[x-2][0]):8.3f}{(atom_coordinate_float[x-2][1]):8.3f}{(atom_coordinate_float[x-2][2]):8.3f}'+ '\n'
                else:
                    line = line
                writer.write(line)