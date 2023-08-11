import functions
import argparse
import os
import Protein
import Ligand 
import AlchemicalFreeEnergy as AFE
import Network
import prepare
import solvate
import logging
logger = logging.getLogger()
logger.setLevel(logging.CRITICAL)


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
    if arguments.working_directory is None:
        arguments.working_directory = os.getcwd()
    elif not os.path.isdir(arguments.working_directory):
        raise argparse.ArgumentTypeError(f"{arguments.working_directory} does not exist")
    
    # type=functions.check_positive(functions.check_int(arguments.min))
    return arguments


def main():
    parser = argparse.ArgumentParser(description="MEZE: MetalloEnZymE FF-builder for alchemistry")

    parser.add_argument("-s",
                        "--step",
                        dest="step",
                        required=True,
                        help="workflow step number: 1: prep, 2: solvate, 3: equilbrate",
                        choices=["1", "2", "3", "4"])

    parser.add_argument("-i",
                        "--ligand-index",
                        dest="idx",
                        type=functions.check_int,
                        help="ligand index for sorting through Network.names",
                        default=0)

    parser.add_argument("-if",
                        "--input-pdb-file",
                        dest="protein",
                        type=functions.file_exists,
                        required=True,
                        help="input pdb file for the metalloenzyme/protein")
    
    parser.add_argument("-g",
                        "--group-name",
                        dest="group_name",
                        help="group name for system, e.g. vim2/kpc2/ndm1/etc",
                        default="meze")

    parser.add_argument("-ppd",
                        "--protein-prep-directory",
                        dest="protein_directory",
                        help="path to protein preparation directory, default: $CWD + /inputs/protein/",
                        default=os.getcwd()+"/inputs/protein/")

    parser.add_argument("-ff",
                        "--force-field",
                        dest="forcefield",
                        help="protein force field, default is ff14SB",
                        default="ff14SB")
    
    parser.add_argument("-w",
                        "--water-model",
                        dest="water_model",
                        help="water model, default is tip3p",
                        default="tip3p")
    
    parser.add_argument("-lpd",
                        "--ligand-prep-directory",
                        dest="ligand_directory",
                        help="path to ligands preparation directory, default: $CWD + /inputs/ligands/",
                        default=os.getcwd()+"/inputs/ligands/")
   
    parser.add_argument("-c",
                        "--ligand-charge",
                        dest="ligand_charge",
                        help="ligand charge",
                        default=0)

    parser.add_argument("-pwd",
                        "--project-working-directory",
                        dest="working_directory",
                        help="working directory for the project",
                        default=os.getcwd())

    parser.add_argument("-e",
                        "--engine",
                        dest="engine",
                        help="MD engine",
                        default="SOMD")

    parser.add_argument("-lf",
                        "--ligand-forcefield",
                        dest="ligand_forcefield",
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
                        default="cubic") 

    parser.add_argument("-st",
                        "--sampling-time",
                        dest="sampling_time",
                        help="sampling time in nanoseconds",
                        default="4") 
    
    # parser.add_argument("-min",
    #                     "--minimisation-steps",
    #                     dest=
    #                     help="number of minimisation steps for equilibration stage",
    #                     default=500)
    
    # parser.add_argument("-snvt",
    #                     "--short-nvt-runtime",
    #                     help="runtime in ps for short NVT equilibration",
    #                     type=functions.check_positive(functions.check_float())) 
    
    # parser.add_argument("-nvt",
    #                     "--nvt-runtime",
    #                     type=functions.check_positive(functions.check_float()),
    #                     help="runtime in ps for NVT equilibration")

    # parser.add_argument("-npt",
    #                     "--npt-runtime",
    #                     type=functions.check_positive(functions.check_float()),
    #                     help="runtime in ps for NPT equilibration",
    #                     default=200)

    arguments = parser.parse_args()
    arguments = clean_arguments(arguments)
    
    protein = Protein.Protein(name=arguments.group_name, 
                              protein_file=arguments.protein, 
                              path=arguments.protein_directory,
                              forcefield=arguments.forcefield,
                              water_model=arguments.water_model)

    network = Network.Network(path=arguments.ligand_directory,
                              forcefield=arguments.ligand_forcefield,
                              charge=arguments.ligand_charge)

    afe = AFE.AlchemicalFreeEnergy(path=arguments.working_directory,
                                   engine=arguments.engine,
                                   sampling_time=arguments.sampling_time,
                                   box_edges=arguments.box_edges,
                                   box_shape=arguments.box_shape)

    if arguments.step == "1":
        prepare.prepare_meze(Protein=protein, Network=network, AFE=afe)
    elif arguments.step == "2":
        solvate.solvate_meze(idx=arguments.idx, Protein=protein, Network=network, AFE=afe)


if __name__ == "__main__":
    main()

