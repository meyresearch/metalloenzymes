import functions
import argparse
import os
import Network
import logging
import time
import equilibrate
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


def main():

    parser = argparse.ArgumentParser(description="MEZE: MetalloEnZymE FF-builder for alchemistry")

    parser.add_argument("protocol_file",
                        help="protocol file containing equilibration options",
                        type=str,
                        default=os.getcwd() + "/afe/protocol.dat")
    
    arguments = parser.parse_args()
    
    protocol = functions.input_to_dict(arguments.protocol_file)

    solvated_network = Network.Network(workdir=protocol["project directory"],
                                       ligand_path=protocol["ligand directory"],
                                       group_name=protocol["group name"],
                                       protein_file=protocol["protein input file"],
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

    equilibrated_network = solvated_network.get_equilibrated()

    equilibrated_network.afe_prep()

    _ = equilibrated_network.write_afe_run_script()
    

if __name__ == "__main__":
    main()

