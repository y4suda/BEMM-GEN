import argparse
import subprocess
import os
import shutil

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
        utils.print_info(f"Amber files are created and saved.")

    return True

def _convert_GROMACS(args: argparse.Namespace):

    parm = parmed.load_file(f"{args.output_prefix}.prmtop", f"{args.output_prefix}.inpcrd")
    parm.save(f"{args.output_prefix}.top", format="gromacs", overwrite=True)
    parm.save(f"{args.output_prefix}.gro", overwrite=True)

    utils.remove_overlap(f"{args.output_prefix}.gro")

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
    print("\n------ amber_simulation.in ------")
    print(f"ntr=1,\nrestraintmask=':{protein_end+1}-{water_start-1} & @C,C1',\nrestraint_wt={restraint_force[0]}")
    print("---------------------------------\n")

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