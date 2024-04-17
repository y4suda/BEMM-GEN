# -*- coding: utf-8 -*-
import argparse
import utils
            

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
    parser_cylinder.add_argument("--length", type=float, default=10.0, help="Length of the cylinder in angstrom. (default: 10.0)")

    # subparser for sphere
    parser_sphere = subparsers.add_parser("sphere", help="Create a sphere model.")

    # add common arguments
    for subparser in [parser_cylinder, parser_sphere]:
        subparser.add_argument("--radius", type=float, default=20.0, help="Radius of the cylinder or sphere in angstrom. (default: 20.0)")
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
        # Error in the options
        if utils.summarize_args(args) is False:
            exit(1)

        if args.command == "cylinder":
            parser_cylinder.print_help()
            ## call cylinder function
        elif args.command == "sphere":
            parser_cylinder.print_help()
            ## call sphere function


    print("Please cite the following paper when you use this software.")
    print("T. Yasuda, R. Morita, Y. Shigeta and R. Harada. (2024) \"Cellular Environment Mimicking Model GENerator: A tool for generating a cellular environment mimicking model.\"")
    