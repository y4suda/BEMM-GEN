import argparse
import utils

def preparation(args: argparse.Namespace):
    
    if validate_args(args) == False:
        raise ValueError("Invalid arguments. Unkown error.")

    if args.proteinpdb is None:
        protein_filename = make_pdb_fromsequence(args.proteinseq, args.proteinSS)
        args.proteinpdb = protein_filename

    major_axis,minor_axis= compute_protein_size(args.proteinpdb)

    radius,length = get_model_size(args, major_axis, minor_axis)
    args.radius = radius
    args.length = length

    return args

def validate_args(args: argparse.Namespace):
    return True

def make_pdb_fromsequence(sequence: str, SS: str):
    protein_filename = "protein.pdb"
    return protein_filename

def compute_protein_size(pdb_filename: str):
    major_axis = None
    minor_axis = None
    return  major_axis,minor_axis

def get_model_size(args: argparse.Namespace, major_axis: float, minor_axis: float):
    radius=None
    length=None
    return radius, length
 

def merge_protein(args: argparse.Namespace):
    #proteinのrorationとtranslationを行う
    return None