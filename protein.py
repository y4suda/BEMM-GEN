import argparse
import subprocess

from sklearn.decomposition import PCA

import utils

AA_DICT = {"A": "ALA",
            "R": "ARG",
            "N": "ASN",
            "D": "ASP",
            "C": "CYS",
            "Q": "GLN",
            "E": "GLU",
            "G": "GLY",
            "H": "HIS",
            "I": "ILE",
            "L": "LEU",
            "K": "LYS",
            "M": "MET",
            "F": "PHE",
            "P": "PRO",
            "S": "SER",
            "T": "THR",
            "W": "TRP",
            "Y": "TYR",
            "V": "VAL"}



def preparation(args: argparse.Namespace):
    
    args = validate_args(args)
    if args == None:
        raise ValueError("Invalid arguments. Unkown error.")

    if args.proteinpdb is None:
        protein_filename = _make_pdb_from_sequence(args.proteinseq, args.proteinSS)
        args.proteinpdb = protein_filename

    major_radius, minor_radius = _compute_protein_size(args.proteinpdb)

    radius, length = _get_model_size(args, major_radius, minor_radius)
    args.radius = radius
    args.length = length

    return args

def validate_args(args: argparse.Namespace):

    if args.proteinpdb is not None and args.proteinseq is not None:
        utils.print_warning("Both protein pdb and sequence are given. Use protein pdb.")

    if args.proteinseq is not None and args.proteinSS is None:
        utils.print_warning("Secondary structure is not given. Use default coil structure.")
        args.proteinSS = "C"*len(args.proteinseq)
    
    if args.proteinSS is not None and len(args.proteinSS) != len(args.proteinseq):
        raise ValueError("Length of the sequence and secondary structure are not same.")
    
    if args.proteinSS is not None:
        for s in args.proteinSS:
            if s not in ["H", "C"]:
                raise ValueError("Secondary structure should be 'H' or 'C'.")

    return args

def _make_pdb_from_sequence(sequence: str, SS: str):
    protein_filename = "protein.pdb"

    utils.print_info(f"Making protein pdb file from the sequence {sequence} and secondary structure {SS}.")
    leap_command = ""
    leap_command += f"source leaprc.protein.ff14SB\n"
    leap_command += f"protein = sequence {{ {' '.join([AA_DICT[aa] for aa in sequence])} }}\n"

    if "H" in SS:
        helix_list = []
        helix_start = None
        helix_end = None
        for i, s in enumerate(SS):
            if s == "H":
                if helix_start is None:
                    helix_start = i
            elif s == "C":
                if helix_start is not None:
                    helix_end = i
                    helix_list.append((helix_start, helix_end))
                    helix_start = None
                    helix_end = None
        if helix_start is not None:
            helix_list.append((helix_start, len(SS)-1))
        
        leap_command += f"impose protein {{ {' '.join([ f'{{ {h[0]} {h[1]} }}' for h in helix_list])} }}   {{ {{ $N $CA  $C $N -40.0 }} {{ $C $N  $CA $C -60.0 }} }}\n"
    leap_command += f"savepdb protein {protein_filename}\n"
    leap_command += f"quit\n"

    with open("leap.in", "w") as f:
        f.write(leap_command)
    
    exitcode = subprocess.run("tleap -f leap.in", shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    if exitcode.returncode != 0:
        raise ValueError("Failed to make protein pdb file in tleap. See leap.log for details.")
    else:
        utils.print_info(f"Protein pdb file is created and saved as {protein_filename}.")

    return protein_filename

def _compute_protein_size(pdb_filename: str):
    major_axis = None
    minor_axis = None
    coords = _load_protein_coordinates(pdb_filename)

    pca = PCA(n_components=3)
    pca.fit(coords)
    coords_pca = pca.transform(coords)

    # compute major and minor axis
    # major_axis = pca.components_[pca.explained_variance_ratio_.argmax()]
    # minor_axis = pca.components_[pca.explained_variance_ratio_.argmin()]

    major_radius = (coords_pca[:,0].max() - coords_pca[:,0].min())/2
    minor_radius = (coords_pca[:,-1].max() - coords_pca[:,-1].min())/2

    utils.print_info(f"Major radius: {major_radius:.2f} Å, Minor radius: {minor_radius:.2f} Å.")
    return  major_radius, minor_radius

def _get_model_size(args: argparse.Namespace, major_radius: float, minor_radius: float):
    radius = None
    length = None

    if args.command == "cylinder":
        # set length
        if args.length is None:
            length = major_radius * 2 + args.padding_length
            utils.print_info(f"Length is not specified. Set to {length:.2f} Å.")
        else:
            length = args.length
        
        if length < major_radius * 2 + args.padding_length:
            utils.print_warning(f"Specified length ({args.length:.2f} Å) is smaller than the protein length ({major_radius:.2f} Å). Protein will be out of the cylinder.")     
        
        # set radius
        if args.radius is None:
            radius = minor_radius + args.padding_radius
            utils.print_info(f"Radius is not specified. Set to {radius:.2f} Å.")
        else:
            radius = args.radius

        if radius < minor_radius + args.padding_radius:
            raise ValueError(f"Specified radius ({args.radius:.2f} Å) is smaller than the protein size ({minor_radius:.2f} Å). Protein will be crushed with the cylinder.")
    elif args.command == "sphere":
        # set radius
        if args.radius is None:
            radius = major_radius + args.padding_radius
            utils.print_info(f"Radius is not specified. Set to {radius:.2f} Å.")
        else:
            radius = args.radius

        if radius < major_radius + args.padding_radius:
            utils.print_warning(f"Specified radius ({args.radius:.2f} Å) is smaller than the protein size ({major_radius:.2f} Å). Protein will be out of the sphere.")
    return radius, length
 

def merge_protein(args: argparse.Namespace):
    #proteinのrorationとtranslationを行う
    return None

def _load_protein_coordinates(pdb_filename: str):
    coords = []
    with open(pdb_filename, "r") as f:
        for line in f:
            if line.startswith("ATOM"):
                coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    return coords