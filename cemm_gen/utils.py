import numpy as np
import sys
import shutil
from scipy.spatial import distance
from logging import basicConfig, getLogger, FileHandler, Formatter, DEBUG, INFO, WARNING, ERROR


def setup_logger(log_level: str):
    basicConfig(level=DEBUG)
    logger = getLogger(__name__)
    logger.propagate = False

    if logger.hasHandlers():
        logger.handlers.clear()

    fh = FileHandler("CEMM-GEN.log")
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

def remove_overlap(filename: str):
    gro_file_name = filename
    shutil.copy(gro_file_name, gro_file_name.replace(".gro","_orig.gro"))

    #座標リストの作成
    atom_coordinate = []
    with open(gro_file_name,"r") as f:
        for x,line in enumerate(f):
            if x > 1 and len(line.split()) != 3:
                atom_coordinate.append(line.split()[-3:])
    atom_coordinate_float = []
    for i in atom_coordinate:
        atom_coordinate_float.append([float(j) for j in i])


    block_size = min(len(atom_coordinate_float), 10000)
    is_overlap =[[True] *((len(atom_coordinate_float)-1)//block_size+1) for i in range((len(atom_coordinate_float)-1)//block_size+1)]

    print_info("Removing atomic clashes...")

    removing_cycle = 10
    for cycle in range(removing_cycle):
        print(f"Starting cycle {cycle}")
        #距離行列の作成
        for block_A in range(0, (len(atom_coordinate_float)-1)//block_size+1):
            for block_B in range(0, (len(atom_coordinate_float)-1)//block_size+1):
                dist = distance.cdist(atom_coordinate_float[block_A*block_size:(block_A+1)*block_size], atom_coordinate_float[block_B*block_size:(block_B+1)*block_size], 'euclidean')
                #オーバーラップの検出
                cutoff = 0.05
                overlap = np.where((dist > 0)&(dist < cutoff))
                overlap_list = list(set(zip(*overlap)))
                overlap_removed = set([str(sorted(i)[0]+block_A*block_size) + "|" + str(sorted(i)[1]+block_B*block_size) for i in overlap_list])
                overlap_removed_list = [[ int(j) for j in (i.split("|"))] for i in overlap_removed]

                if len(overlap_removed_list) != 0:
                    is_overlap[block_A][block_B] = True
                else:
                    is_overlap[block_A][block_B] = False

                for j in overlap_removed_list:
                    atom_a = np.array(atom_coordinate_float[j[0]])
                    atom_b = np.array(atom_coordinate_float[j[1]])
                    vector = atom_a - atom_b
                    distance_a_b = np.linalg.norm(vector)
                    if distance_a_b != 0:
                        atom_b_n = atom_b - vector/distance_a_b * (cutoff/2)
                    else:
                        atom_b_n = atom_b + np.random.rand(3) * (cutoff/2)
                    atom_coordinate_float[j[1]] = atom_b_n
        
        if not any(np.array(is_overlap).flatten()):
            print_info("Complete removing atomic crashes. There is no overlap.") 
            break
        if cycle == removing_cycle-1:
            print_warning(f"There are still overlaps after {removing_cycle} cycles. Please check the output file carefully.")

    with open(gro_file_name, "w") as writer:
        with open(gro_file_name.replace(".gro","_orig.gro"), "r") as f:
            for x, line in enumerate(f):
                if x > 1 and len(line.split()) != 3 and "WAT" not in line and "SOL" not in line:
                    line = line[:20] +f'{(atom_coordinate_float[x-2][0]):8.3f}{(atom_coordinate_float[x-2][1]):8.3f}{(atom_coordinate_float[x-2][2]):8.3f}'+ '\n'
                else:
                    line = line
                writer.write(line)