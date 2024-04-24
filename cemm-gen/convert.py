import argparse
import subprocess
import os
import shutil

import parmed

import utils


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
    
    # for resname in args.resnames.split(":"):
    #     leap_command += f"loadamberparams {resname}.frcmod\n"
    #     leap_command += f"loadoff {resname}.lib\n"

    # if args.command == "cylinder":
    #     leap_command += f"model = loadpdb cylinder.pdb\n"
    # elif args.command == "sphere":
    #     leap_command += f"model = loadpdb sphere.pdb\n"

    if args.proteinpdb is not None:
        leap_command += f"protein = loadpdb {args.proteinpdb}\n"
        leap_command += f"system = protein\n"
        # leap_command += f"system = combine {{ model protein }}\n"
    else:
        leap_command += f"system = model\n"

    leap_command += f"solvateBox system TIP3PBOX 10 \n"
    leap_command += f"addions system Na+ 0\n"
    leap_command += f"addions system Cl- 0\n"
    leap_command += f"charge system\n"

    leap_command += f"saveamberparm system {args.output_prefix}.prmtop {args.output_prefix}.inpcrd\n"
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

    utils.print_info(f"GROMACS files are created and saved.")
    return True

def _restraint_amber(args: argparse.Namespace):
    utils.print_info("Use the following options in .in files for the Amber simulation.")
    model_resn = 0
    with open(args.proteinpdb, "r") as f:
        for line in f:
            if line.startswith("ATOM") and (" C " in line[12:15] or "C1" in line[12:15]):
                model_resn += 1

    restraint_force = [1000, 1000, 1000]
    print("\n------ amber_simulation.in ------")
    print(f"ntr=1,\nrestraintmask=':1-{model_resn} & C1',\nrestraint_wt={restraint_force[0]}")
    print("---------------------------------\n")

    return True

def _restraint_GROMACS(args: argparse.Namespace):
    
    # get residue id of CA atoms
    restraint_resid = []
    with open(f"{args.output_prefix}.gro", "r") as f:
        for line in f:
            if " C " in line[10:15] or "C1" in line[10:15]:
                restraint_resid.append(int(line[15:20]))

    # write restraint.itp
    restraint_force = [1000, 1000, 1000]
    with open("restraint.itp", "w") as f:
        f.write("[ position_restraints ]\n")
        for resn in restraint_resid:
            f.write(f"{resn: >5}    1    {restraint_force[0]} {restraint_force[1]} {restraint_force[2]}\n")

    # add restraint.itp to the top file
    shutil.move(f"{args.output_prefix}.top", f"{args.output_prefix}_norestraint.top")
    moleculetype_count = 0
    with open(f"{args.output_prefix}_norestraint.top", "r") as f:
        with open(f"{args.output_prefix}.top", "w") as g:
            for line in f:
                if line.startswith("[ moleculetype ]"):
                    moleculetype_count += 1
                    if moleculetype_count == 2:
                        g.write("#include \"restraint.itp\"\n\n")
                g.write(line)
    os.remove(f"{args.output_prefix}_norestraint.top")
    utils.print_info("Position restraints were added in topology file for the GROMACS simulation.")

    return True