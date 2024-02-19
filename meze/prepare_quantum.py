
import argparse 
import os
import functions
from qmMeze import Meze



# def get_default_config(config_file):



def main():

    parser = argparse.ArgumentParser(description="solvation for meze workflow")

    parser.add_argument("ligand_name",
                        help="ligand name",
                        type=str)
    
    parser.add_argument("protocol_file",
                        help="protocol file containing equilibration options",
                        type=str,
                        default=os.getcwd() + "/qmmm_input_files/protocol.dat")

    arguments = parser.parse_args()
    protocol = functions.input_to_dict(arguments.protocol_file)

    meze = Meze(workdir=protocol["project directory"],
                is_qm=True,
                # ligand_path=protocol["ligand directory"],
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
                pressure=protocol["pressure"],
                output=protocol["outputs"])

    solvated_meze = meze.set_universe(file_name=meze.protein_path + "bound_" + arguments.ligand_name + "_solvated")
    # write scripts
    solvated_meze.model_0(ligand_name=arguments.ligand_name)
    

if "__name__" == main():
    main()



