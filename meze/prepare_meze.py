import functions
import argparse
import os

def clean_arguments(arguments):
    """
    Check arguments and clean them

    Parameters:
    -----------
    arguments: Namespace
        command-line arguments
    Return:
    -------
    cleaned Namespace of arguments
    
    """
    if arguments.workdir is None:
        arguments.workdir = os.getcwd()
    elif not os.path.isdir(arguments.workdir):
        raise argparse.ArgumentTypeError(f"{arguments.workdir} does not exist")
    return arguments


def main():
    parser = argparse.ArgumentParser(description="MEZE: MetalloEnZymE FF-builder for alchemistry")
   
    parser.add_argument("-i",
                        "--input-pdb",
                        dest="protein",
                        required=True, 
                        type=functions.file_exists,
                        help="input pdb file for the metalloenzyme/protein")
    
    parser.add_argument("-g",
                        "--group-name",
                        dest="group_name",
                        help="group name for system, e.g. vim2/kpc2/ndm1/etc",
                        default="meze")

    parser.add_argument("-wd",
                        "--working_directory",
                        dest="workdir",
                        help="project working directory",
                        default=os.getcwd())
    
    parser.add_argument("-ff",
                        "--force-field",
                        dest="forcefield",
                        help="protein force field",
                        default="ff14SB")
    
    parser.add_argument("-w",
                        "--water-model",
                        dest="water_model",
                        help="water model",
                        default="tip3p")
    
    parser.add_argument("-e",
                        "--engine",
                        dest="engine",
                        help="MD engine",
                        default="SOMD")

    parser.add_argument("-lf",
                        "--ligand-forcefield",
                        dest="ligand_ff",
                        help="ligand force field",
                        default="gaff2")
    
    parser.add_argument("-be",
                        "--box-edges",
                        dest="box_edges",
                        help="box edges in angstroms",
                        default="20") 
    
    parser.add_argument("-bs",
                        "--box-shape",
                        dest="box_shape",
                        help="box shape",
                        default="orthorhombic") 

    parser.add_argument("-st",
                        "--sampling-time",
                        dest="sampling_time",
                        help="sampling time in nanoseconds",
                        default="4") 

    arguments = parser.parse_args()
    arguments = clean_arguments(arguments)
    
    # this should be somehow an object where all the arguments are attributes 
    # then i can access them from every function in an easy way
    functions.prepare_system(arguments.group_name,
                             arguments.workdir, 
                             arguments.protein, 
                             arguments.forcefield,
                             arguments.water_model,
                             arguments.engine,
                             arguments.ligand_ff,
                             arguments.box_edges,
                             arguments.box_shape,
                             arguments.sampling_time)

if __name__ == "__main__":
    main()

