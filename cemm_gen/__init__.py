# -*- coding: utf-8 -*-
import os
import sys
import argparse
from importlib.metadata import version, PackageNotFoundError
from . import utils
from . import protein
from . import model
from . import convert
from . import make_param

def main():
    
    try:
        __version__ = version("your-package-name")
    except PackageNotFoundError:
        __version__ = "0.0.0"

    logger = utils.setup_logger("INFO")
    logger.info(f" ")
    logger.info(f"Start CEMM-GEN. Version: {__version__}")
    logger.info("Command line arguments: " + " ".join(sys.argv))

    print(f"""
          
 ██████╗███████╗███╗   ███╗███╗   ███╗       ██████╗ ███████╗███╗   ██╗
██╔════╝██╔════╝████╗ ████║████╗ ████║      ██╔════╝ ██╔════╝████╗  ██║
██║     █████╗  ██╔████╔██║██╔████╔██║█████╗██║  ███╗█████╗  ██╔██╗ ██║
██║     ██╔══╝  ██║╚██╔╝██║██║╚██╔╝██║╚════╝██║   ██║██╔══╝  ██║╚██╗██║
╚██████╗███████╗██║ ╚═╝ ██║██║ ╚═╝ ██║      ╚██████╔╝███████╗██║ ╚████║
 ╚═════╝╚══════╝╚═╝     ╚═╝╚═╝     ╚═╝       ╚═════╝ ╚══════╝╚═╝  ╚═══╝
          
          Cellular Environment Mimicking Model GENerator  {__version__}

                T. Yasuda, R. Morita, Y. Shigeta and R. Harada. (2024)
          
          """)

    # main parser
    parser= argparse.ArgumentParser(description="Make a cellular environment mimicking model.", add_help=False)
    subparsers = parser.add_subparsers(dest="command")

    # subparser for cylinder
    parser_cylinder = subparsers.add_parser("cylinder", help="Create a cylinder model.")
    parser_cylinder.add_argument("--length", type=float, help="Length of the cylinder in angstrom.")
    parser_cylinder.add_argument("--padding-length", type=float, default=10.0, help="Padding length in angstrom. (default: 10.0)")
    # subparser for sphere
    parser_sphere = subparsers.add_parser("sphere", help="Create a sphere model.")

    # add common arguments
    for subparser in [parser_cylinder, parser_sphere]:
        subparser.add_argument("--radius", type=float, help="Radius of the cylinder or sphere in angstrom.")
        subparser.add_argument("--padding-radius", type=float, default=10.0, help="Padding radius in angstrom. (default: 10.0)")

        subparser.add_argument("--proteinpdb", type=str, default=None, help="PDB file of the protein to be placed in the model.")
        subparser.add_argument("--proteinseq", type=str, default=None, help="Amino acid sequence of the protein to be placed in the model.")
        subparser.add_argument("--proteinSS", type=str, default=None, help="Secondary structure of the protein to be placed in the model.")

        subparser.add_argument("--forcefield-protein", type=str, default="ff14SB", choices=["ff14SB", "ff99SB"], help="Force field for the protein. (default: ff14SB)")
        subparser.add_argument("--forcefield-model", type=str, default="gaff2", choices=["gaff2", "gaff"], help="Force field for the residues. (default: gaff2)")

        subparser.add_argument("--min-distance", type=float, default=4.0, help="Minimum distance between two residues in angstrom. (default: 4.0)")
        subparser.add_argument("--resnames", type=str, default="BOC", help="Residue names to use connected by colon. Ex. \"BOC:BOH\" (default: BOC)")
        subparser.add_argument("--composition", type=str, default="1", help="Composition of the residues in ratio connected by colon, raito or percentage is acceptable. Ex. \"1:1\" (default: 1)")
        
        subparser.add_argument("--outward", action="store_true", help="Place residues outward from the center. (dafault: False)")
        subparser.add_argument("--output-prefix", type=str, default="out", help="Prefix for output file name. (default: inward)")

    parser_makeparam = subparsers.add_parser("makeparam", help="Make parameter files for residues in model.")
    parser_makeparam.add_argument("--smiles", type=str, help="SMILES file.")
    parser_makeparam.add_argument("--resname", required=True, type=str, help="Residue name.")
    parser_makeparam.add_argument("--description", required=True, type=str, default="Description", help="Description of the residue.")
    parser_makeparam.add_argument("--overwrite", action=argparse.BooleanOptionalAction, default=False, help="Overwrite the existing parameter files. (defailt: False)")
    parser_makeparam.add_argument("--method", type=str, default="HF", help="Method for RESP calculation.")
    parser_makeparam.add_argument("--basisSet", type=str, default="6-31G*", help="Basis set for RESP calculation.")
    parser_makeparam.add_argument("--method-opt", type=str, default="B3LYP", help="Method for optimization.")
    parser_makeparam.add_argument("--basisSet-opt", type=str, default="6-31G*", help="Basis set for optimization.")
    parser_makeparam.add_argument("--num-thread", type=int, default=8, help="Number of threads.")
    parser_makeparam.add_argument("--memory-sizeGB", type=str, default="8", help="Memory size.")
    parser_makeparam.add_argument("--neutralize", action=argparse.BooleanOptionalAction, default=True, help="Neutralize the molecule. (defailt: False)")
    parser_makeparam.add_argument("--singlePoint", action=argparse.BooleanOptionalAction, default=False, help="Single point calculation. (defailt: False)")
    parser_makeparam.add_argument("--netcharge", type=int, default=0, help="Net charge of the molecule.")
    parser_makeparam.add_argument("--multiplicity", type=int, default=1, help="Multiplicity of the molecule.")
    parser_makeparam.add_argument("--max-iter", type=int, default=500, help="Maximum number of iterations for optimization.")
    parser_makeparam.add_argument("--peptideseq", type=str, default=None, help="Peptide sequence.")

    parser_listparam = subparsers.add_parser("listparam", help="List parameter files for residues in model.")
    parser_listparam.add_argument("--dump", action=argparse.BooleanOptionalAction, default=False, help="Dump the parameter files.")
    args = parser.parse_args()


    # No options are given, print help
    if args.command is None:
        
        parser.print_help()
    elif args.command == "makeparam":
        utils.print_info("Make parameter files for residues in model...")
        make_param.make_param(args)
    elif args.command == "listparam":
        utils.print_info("List parameter files for residues in model...")
        make_param.list_param(args)
    else:

        if args.proteinpdb is not None or args.proteinseq is not None or args.proteinSS is not None:
            utils.print_info("Prepare protein...")
            protein.preparation(args)

        if args.command == "cylinder":
            ## call cylinder function
            utils.print_info("Create a cylinder model...")
            model.build_cylinder(args)
        elif args.command == "sphere":
            ## call sphere function
            utils.print_info("Create a sphere model...")
            model.build_sphere(args)
        
        if args.proteinpdb is not None:

            protein.centering_rotation_protein(args)

        utils.summarize_args(args)
        utils.print_info("Convert the model to MD input files...")
        convert.make_MD_input(args)

    print("\n")
    print("Please cite the following paper when you use this software.")
    print("T. Yasuda, R. Morita, Y. Shigeta and R. Harada. (2024) \"Cellular Environment Mimicking Model GENerator: A tool for generating a cellular environment mimicking model.\"")
    print("")
    
    logger.info("End CEMM-GEN.")