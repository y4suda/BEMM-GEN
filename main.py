# -*- coding: utf-8 -*-
import argparse
import utils
import protein
import model
import convert

if __name__ == "__main__":

    print("""
          
 ██████╗███████╗███╗   ███╗███╗   ███╗       ██████╗ ███████╗███╗   ██╗
██╔════╝██╔════╝████╗ ████║████╗ ████║      ██╔════╝ ██╔════╝████╗  ██║
██║     █████╗  ██╔████╔██║██╔████╔██║█████╗██║  ███╗█████╗  ██╔██╗ ██║
██║     ██╔══╝  ██║╚██╔╝██║██║╚██╔╝██║╚════╝██║   ██║██╔══╝  ██║╚██╗██║
╚██████╗███████╗██║ ╚═╝ ██║██║ ╚═╝ ██║      ╚██████╔╝███████╗██║ ╚████║
 ╚═════╝╚══════╝╚═╝     ╚═╝╚═╝     ╚═╝       ╚═════╝ ╚══════╝╚═╝  ╚═══╝
          
          Cellular Environment Mimicking Model GENerator  ver 0.01

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

        subparser.add_argument("--forcefield-protein", type=str, default="ff14SB", help="Force field for the protein. (default: ff14SB)")
        subparser.add_argument("--forcefield-model", type=str, default="gaff2", help="Force field for the residues. (default: gaff2)")

        subparser.add_argument("--min-distance", type=float, default=4.0, help="Minimum distance between two residues in angstrom. (default: 4.0)")
        subparser.add_argument("--resnames", type=str, default="BOC", help="Residue names to use connected by colon. Ex. \"BOC:BOH\" (default: BOC)")
        subparser.add_argument("--composition", type=str, default="1", help="Composition of the residues in ratio connected by colon, raito or percentage is acceptable. Ex. \"1:1\" (default: 1)")
        
        subparser.add_argument("--outward", action="store_true", help="Place residues outward from the center. (dafault: False)")
        subparser.add_argument("--output-prefix", type=str, default="out", help="Prefix for output file name. (default: inward)")

    # parse arguments
    args = parser.parse_args()


    # No options are given, print help
    if args.command is None:
        
        parser.print_help()
    else:

        if args.proteinpdb is not None or args.proteinseq is not None or args.proteinSS is not None:
            protein.preparation(args)

        if args.command == "cylinder":
            ## call cylinder function
            model.build_cylinder(args)
        elif args.command == "sphere":
            ## call sphere function
            model.build_sphere(args)
        
        if args.proteinpdb is not None:
            protein.centering_rotation_protein(args)

        utils.summarize_args(args)
        convert.make_MD_input(args)

    print("\n")
    print("Please cite the following paper when you use this software.")
    print("T. Yasuda, R. Morita, Y. Shigeta and R. Harada. (2024) \"Cellular Environment Mimicking Model GENerator: A tool for generating a cellular environment mimicking model.\"")
    print("")
    