import sys
import shutil
from logging import basicConfig, getLogger, FileHandler, Formatter, DEBUG, INFO, WARNING, ERROR


def setup_logger(log_level: str):
    basicConfig(level=DEBUG)
    logger = getLogger(__name__)
    logger.propagate = False

    if logger.hasHandlers():
        logger.handlers.clear()

    fh = FileHandler("BEMM-GEN.log")
    fh.setLevel(DEBUG)
    formatter = Formatter("%(asctime)s - %(levelname)s - %(message)s")
    fh.setFormatter(formatter)
    logger.addHandler(fh)
    return logger

logger = setup_logger("INFO")


def red_text(text: str):
    """
    Return red colored text.
    """
    return f"\033[91m{text}\033[0m"

def orange_text(text: str):
    """
    Return orange colored text.
    """
    return f"\033[93m{text}\033[0m"

def green_text(text: str):
    """
    Return green colored text.
    """
    return f"\033[92m{text}\033[0m"

def print_error(text: str):
    """
    Print error message.
    """
    print(f"{red_text('Error')}: {text}")
    logger.error(text)
    return

def print_warning(text: str):
    """
    Print warning message.
    """
    print(f"{orange_text('Warning')}: {text}")
    logger.warning(text)
    return

def print_info(text: str):
    """
    Print information message.
    """
    print(f"{green_text('INFO')}: {text}")
    logger.info(text)
    return True

def print_debug(text: str):
    """
    Print debug message.
    """
    logger.debug(text)
    return True

def summarize_args(args):
    """
    Summarize arguments.
    """

    print("\nSummary of user options\n")
    if "command" in args:
        print(f"Model type:\t\t{args.command}")
        logger.info(f"Model type:\t\t{args.command}")
    if "length" in args and args.length is not None:
        print(f"Length:\t\t{args.length:.2f} Å")
        logger.info(f"Length:\t\t{args.length:.2f} Å")
    if "radius" in args:
        print(f"Radius:\t\t{args.radius:.2f} Å")
        logger.info(f"Radius:\t\t{args.radius:.2f} Å")
    if "min_distance" in args:
        print(f"Min. distance:\t{args.min_distance} Å")
        logger.info(f"Min. distance:\t{args.min_distance} Å")
    
    #この処理はmodel.pyに移動する
    if "resnames" in args and "composition" in args:
        num_resnames = len(args.resnames.split(":"))
        num_composition = len(args.composition.split(":"))
        if  num_resnames != num_composition:
            print_error(f"Number of resnames ({num_resnames}, {args.resnames}) and composition ({num_composition}, {args.composition}) are not same.")
            return False
        else:
            print(f"Residue names:\t{args.resnames.split(':')}")
            print(f"Composition:\t{list(map(float, args.composition.split(':')))}")

            
    if "outward" in args:
        print(f"Outward:\t{args.outward}")
    if "output_prefix" in args:
        print(f"Output prefix:\t{args.output_prefix}")

    print("\nInputs are validated.\n")
    return True