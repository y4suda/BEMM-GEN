import os
import sys
import glob
import psi4
import resp
import rdkit
import shutil
import argparse
import subprocess


from openbabel import openbabel as ob
from rdkit import Chem
from rdkit.Chem import AllChem,Draw

import utils

# This code is partialy based on the following repository:
# https://github.com/pablo-arantes/making-it-rain
# Copyright (c) 2021 Pablo Ricardo Arantes

AA_list =  ["ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL", "HIE", "HID", "HIP", "CYM", "CYX"]


def _neutralize_atoms(mol):
    pattern = Chem.MolFromSmarts("[+1!h0!$([*]~[-1,-2,-3,-4]),-1!$([*]~[+1,+2,+3,+4])]")
    at_matches = mol.GetSubstructMatches(pattern)
    at_matches_list = [y[0] for y in at_matches]
    if len(at_matches_list) > 0:
        for at_idx in at_matches_list:
            atom = mol.GetAtomWithIdx(at_idx)
            chg = atom.GetFormalCharge()
            hcount = atom.GetTotalNumHs()
            atom.SetFormalCharge(0)
            atom.SetNumExplicitHs(hcount - chg)
            atom.UpdatePropertyCache()
    return mol

def _cleanUp(psi4out_xyz):
    deleteTheseFiles = ["1_default_grid.dat","1_default_grid_esp.dat","grid.dat","timer.dat"]
    deleteTheseFiles.append(psi4out_xyz)
    for fileName in deleteTheseFiles:
        if os.path.exists(fileName):
            os.remove(fileName)

def _get_xyz_coords(mol):
    if not mol is None:
        num_atoms = mol.GetNumAtoms()
        xyz_string=""
        for counter in range(num_atoms):
            pos=mol.GetConformer().GetAtomPosition(counter)
            xyz_string = xyz_string + ("%s %12.6f %12.6f %12.6f\n" % (mol.GetAtomWithIdx(counter).GetSymbol(), pos.x, pos.y, pos.z) )
    return xyz_string


def _calcRESPCharges(mol, basisSet, method, gridPsi4 = 16):
    options = {"BASIS_ESP": basisSet,
               "METHOD_ESP": method,
               "RESP_A": 0.0005,
               "RESP_B": 0.1,
               "VDW_SCALE_FACTORS":[1.4, 1.6, 1.8, 2.0],
               "VDW_POINT_DENSITY":int(gridPsi4)
    }

    # resp_charges = resp.resp([mol], [options])[0][1]
    resp_charges = resp.resp([mol], options)

    return resp_charges


def calc_RESP_charges(mol_smiles, resname, method="HF", basisSet="6-31G*", method_opt="B3LYP", basisSet_opt="6-31G*", neutralize=True, singlePoint=True, num_thread=8, memory_sizeGB="8 GB", path="./", netcharge=0, multiplicity=1):

    psi4.set_num_threads(num_thread)
    psi4.set_memory(memory_sizeGB)
    psi4.set_output_file("psi4.log")
    psi4.set_options({"g_convergence": "gau"})
    obConversion = ob.OBConversion()
    obConversion.SetInAndOutFormats("xyz", "mol2")

    SMILESasInput = True

    mol_noH = Chem.MolFromSmiles(mol_smiles)
    mol = Chem.AddHs(mol_noH, addCoords=True)

    Draw.MolToFile(mol_noH, size=(600, 600), filename=f"{path}/smiles.png")

    try:
        # MMFF optimization
        # AllChem.MMFFOptimizeMolecule(mol)

        # ETKDGv3 optimization
        p = AllChem.ETKDGv3()
        AllChem.EmbedMolecule(mol, p)
    except:
        print("ETKDGv3 failed, trying UFF optimization")
        AllChem.UFFOptimizeMolecule(mol)

    Draw.MolToFile(mol, size=(600, 600), filename=f"{path}/smiles_ETKDGv3.png")
    
    if not mol is None:
        # molId = mol.GetProp("_Name")
        molId = resname
        print("Trying:", molId)

        if neutralize:
            mol = _neutralize_atoms(mol)
            mol = Chem.AddHs(mol)

        charge_string = f"{netcharge} {multiplicity}\n" 
        xyz_string = _get_xyz_coords(mol)
        psi_mol = psi4.geometry(charge_string+xyz_string)

        ### single point calculation
        outfile_mol2 = resname+".mol2"

        if singlePoint:
            print("Running singlepoint...")
            resp_charges = _calcRESPCharges(psi_mol, basisSet, method, gridPsi4 = 1)

        else:
            print("Running geometry optimization...")
            methodNbasisSet = method_opt+"/"+basisSet_opt
            psi4.optimize(methodNbasisSet, molecule=psi_mol)
            resp_charges = _calcRESPCharges(psi_mol, basisSet, method, gridPsi4 = 1)

        ### save coords to xyz file
        psi4out_xyz = molId + ".xyz"
        psi_mol.save_xyz_file(psi4out_xyz,1)


        ### read xyz file and write as mol2
        ob_mol = ob.OBMol()
        obConversion.ReadFile(ob_mol, psi4out_xyz)

        ### write as mol2
        outfile_mol2 = path+"/"+molId+"_partialChgs.mol2"
        obConversion.WriteFile(ob_mol, outfile_mol2)

        ### set new partial charges
        count = 0
        newChg_temp = resp_charges[1]
        # print("RESP Charges: ", newChg_temp)
        for atom in ob.OBMolAtomIter(ob_mol):
            newChg = newChg_temp[count]
            atom.SetPartialCharge(newChg)
            count += 1

        ### write as mol2
        outfile_mol2 = path+"/"+molId+".mol2"
        outfile_pdb = path+"/"+molId+".pdb"
        utils.print_info("Finished. Saved compound with partial charges as mol2 file: %s" % outfile_mol2)
        obConversion.WriteFile(ob_mol, outfile_mol2)
        ## clean up
        _cleanUp(psi4out_xyz)

    #draw_with_charges

    for at, i in zip(mol.GetAtoms(), newChg_temp):
        lbl = "%.2f"%(i)
        at.SetProp("atomNote",lbl)
        Draw.MolToFile(mol, size=(600, 600), filename=f"{path}/smiles_charges.png")

def _round_mol2_charge(mol2_file):
    # Mol2ファイルの読み込み
    with open(mol2_file, "r") as f:
        lines = f.readlines()

    # 電荷の丸め込みと合計の算出
    sum_charge = 0.0
    charges = []
    atom_types = []
    for i, line in enumerate(lines):
        if line.startswith("@<TRIPOS>ATOM"):
            for l in lines[i+1:]:
                if l.startswith("@<TRIPOS>"):
                    break
                charges.append(round(float(l.split()[8]), 5))
                atom_types.append(l.split()[1])

    print(f"Rounded charges: {charges}")
    sum_charge = round(sum(charges), 5)

    # 電荷の整数部分の算出

    # 電荷の差分の算出と適用
    diff_charge = round(round(sum_charge,0) - sum_charge, 5)
    mod_charge_atom_index = -1
    for index, at in enumerate(atom_types):
        if "H" in at:
            continue
        else:
            mod_charge_atom_index = index
            
    # 適当な原子に電荷を加える
    charges[mod_charge_atom_index] += diff_charge
    print(f"Modified charges: {charges}")

    # 電荷の合計の再算出
    sum_charge = sum(charges)

    # 調整した結果の出力
    is_atom = False
    atom_index = 0
    with open(mol2_file.replace(".mol2", "_round.mol2"), "w") as f:
        for i, line in enumerate(lines):
            if line.startswith("@<TRIPOS>ATOM"):
                is_atom = True
                f.write(line)
                continue
            if line.startswith("@<TRIPOS>"):
                is_atom = False
                f.write(line)
                continue
            if is_atom == True:
                parts = line.split()
                parts[8] = f"{round(charges[atom_index], 5): <2.5f}"
                atom_index += 1
                f.write("\t".join(parts) + "\n")
                continue
            f.write(line)



def make_param(args : argparse.Namespace):

    smiles = args.smiles

    if len(args.resname) > 3:
        utils.print_error("Residue name should be less than 3 characters.")
        sys.exit(1)

    if args.resname in AA_list:
        utils.print_error("Residue name is already used in AMBER force field.")
        sys.exit(1)

    if not args.smiles.startswith("CC"):
        utils.print_error("SMILES should start with 'CC'.")
        sys.exit(1)


    res_dir = f"./FF_PARAM/FF_{args.resname}"
    os.makedirs(f"{res_dir}", exist_ok=True)

    with open(f"{res_dir}/description.txt", "w") as f:
        f.write(args.description)

    calc_RESP_charges(smiles, resname=args.resname , path=f"{res_dir}", neutralize=args.neutralize, singlePoint=args.singlePoint, num_thread=args.num_thread, memory_sizeGB=f"{args.memory_sizeGB} GB",
                        method=args.method, basisSet=args.basisSet, method_opt=args.method_opt, basisSet_opt=args.basisSet_opt, netcharge=args.netcharge, multiplicity=args.multiplicity)

    
    if not os.path.exists(f"{res_dir}/{args.resname}.mol2"):
        utils.print_error("Failed to calculate RESP charges.")
        sys.exit(1)

    utils.print_info(f"RESP charges are calculated and saved as {res_dir}/{args.resname}.mol2.")

    os.remove(f"./psi4.log")

    subprocess.run(f"antechamber -i {res_dir}/{args.resname}.mol2 -fi mol2 -o {res_dir}/{args.resname}_orig.mol2 -fo mol2 -at gaff2 -rn {args.resname} -nc {args.netcharge} -m {args.multiplicity} -pf y", shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    _round_mol2_charge(f"{res_dir}/{args.resname}_orig.mol2")
    shutil.copyfile(f"{res_dir}/{args.resname}_orig_round.mol2", f"{res_dir}/{args.resname}.mol2")
    subprocess.run(f"parmchk2 -i {res_dir}/{args.resname}.mol2 -f mol2 -o {res_dir}/{args.resname}.frcmod", shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    # os.remove("./timer.dat")
    os.remove("./results.out")

    leap_command = ""
    leap_command += f"source leaprc.protein.ff14SB\n"
    leap_command += f"source leaprc.gaff2\n"
    leap_command += f"{args.resname} = loadmol2 {res_dir}/{args.resname}.mol2\n"
    leap_command += f"loadamberparams {res_dir}/{args.resname}.frcmod\n"
    leap_command += f"saveoff {args.resname} {res_dir}/{args.resname}.lib\n"
    leap_command += f"savepdb {args.resname} {res_dir}/{args.resname}.pdb\n"
    leap_command += f"quit"

    with open(f"{res_dir}/leap.in", "w") as f:
        f.write(leap_command)

    exitcode = subprocess.run(f"tleap -f {res_dir}/leap.in", shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    if exitcode.returncode != 0:
        utils.print_error(f"Failed to generate force field parameter for {args.resname} in tleap. See leap.log for details.")
        exit(1)
    else:
        utils.print_info(f"Force field parameter files are created and saved.")
        os.remove(f"./leap.log")

    return True


def list_param(args : argparse.Namespace):
    print("")
    print(f"name: description")
    print(f"{'-'*60}")
    for d in glob.glob("./FF_PARAM/FF_*"):
        resname = d.split("/")[-1].split("_")[-1]
        with open(f"{d}/description.txt", "r") as f:
            description = f.read()
        print(f"{resname: <4}: {description}")

    print("")
    return True