import functions
import argparse
import os
import sofra
import logging
import meze.meze as meze
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
    check_integers = lambda arg: functions.check_positive(functions.check_int(arg))
    check_integers(arguments.min_steps)
    check_floats = lambda arg: functions.check_positive(functions.check_float(arg))
    check_floats(arguments.short_nvt)
    check_floats(arguments.nvt)
    check_floats(arguments.npt)

    return arguments


def character(string):
    if len(string) > 1:
        message = f"Input must be a single character, not {string}"
        raise argparse.ArgumentTypeError(message)
    if string == "_":
        message = "The character _ is not an allowed separator"
        raise argparse.ArgumentTypeError(message)


def main():

    parser = argparse.ArgumentParser(description="MEZE: MetalloEnZymE FF-builder for alchemistry")

    parser.add_argument("protocol_file",
                        help="protocol file containing equilibration options",
                        type=str,
                        default=os.getcwd() + "/afe/protocol.dat")
    
    parser.add_argument("transformation",
                        help="the pair of ligands undergoing AFE transformation, e.g. ligand_1~ligand_2",
                        type=str)
    
    parser.add_argument("-s",
                        "--separator",
                        help="character separating the two ligand names",
                        default="~",
                        type=character)
    
    arguments = parser.parse_args()
    
    protocol = functions.input_to_dict(arguments.protocol_file)

    keys = list(protocol.keys())
    if "metal" in keys:
        metal = True
    else:
        metal = False
    
    if metal:
        meze = meze.Meze(prepared=True,
                            protein_file=protocol["protein input file"],
                            cut_off=protocol["cutoff"],
                            force_constant_0=protocol["force constant"],
                            workdir=protocol["project directory"],
                            ligand_path=protocol["ligand directory"],
                            ligand_charge=protocol["ligand charge"],
                            ligand_ff=protocol["ligand forcefield"],
                            group_name=protocol["group name"],        
                            protein_path=protocol["protein directory"],
                            water_model=protocol["water model"],
                            protein_ff=protocol["protein forcefield"],
                            engine=protocol["engine"],
                            sampling_time=protocol["engine"],
                            box_edges=protocol["box edges"],
                            box_shape=protocol["box shape"],
                            min_steps=protocol["minimisation steps"],
                            short_nvt=protocol["short nvt"],
                            nvt=protocol["nvt"],
                            npt=protocol["npt"],
                            min_dt=protocol["minimisation stepsize"],
                            min_tol=protocol["minimisation tolerance"],
                            repeats=protocol["repeats"],
                            temperature=protocol["temperature"],
                            pressure=protocol["pressure"])
        
        meze.write_restraints_file_0(engine="somd")

    elif not metal:

        meze = sofra.Sofra(prepared=True,
                               equilibration_path=protocol["equilibration directory"],
                               outputs=protocol["outputs"],
                               workdir=protocol["project directory"],
                               ligand_path=protocol["ligand directory"],
                               group_name=protocol["group name"],
                               protein_file=protocol["prepared protein file"],
                               protein_path=protocol["protein directory"],
                               water_model=protocol["water model"],
                               ligand_ff=protocol["ligand forcefield"],
                               protein_ff=protocol["protein forcefield"],
                               ligand_charge=protocol["ligand charge"],
                               engine=protocol["engine"],
                               sampling_time=protocol["sampling time"],
                               box_edges=protocol["box edges"],
                               box_shape=protocol["box shape"],
                               min_steps=protocol["minimisation steps"],
                               short_nvt=protocol["short nvt"],
                               nvt=protocol["nvt"],
                               npt=protocol["npt"],
                               min_dt=protocol["minimisation stepsize"],
                               min_tol=protocol["minimisation tolerance"],
                               repeats=protocol["repeats"],
                               temperature=protocol["temperature"],
                               pressure=protocol["pressure"])
          
    ligand_a, ligand_b = functions.separate(arguments.transformation)

    equilibrated_network = meze.get_equilibrated(ligand_a, ligand_b)
    equilibrated_network.prepare_afe(ligand_a, ligand_b) 

if __name__ == "__main__":
    main()

