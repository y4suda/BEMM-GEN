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
    _convert_GROMACS(args)

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
        raise ValueError("Failed to convert pdb file to amber format in tleap. See leap.log for details.")
    else:
        utils.print_info(f"Amber files are created and saved.")

    return True

def _convert_GROMACS(args: argparse.Namespace):

    parm = parmed.load_file(f"{args.output_prefix}.prmtop", f"{args.output_prefix}.inpcrd")
    parm.save(f"{args.output_prefix}.top", format="gromacs", overwrite=True)
    parm.save(f"{args.output_prefix}.gro", overwrite=True)

    utils.print_info(f"GROMACS files are created and saved.")
    return True
