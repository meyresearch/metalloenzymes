import argparse
import functions
import os
import Network
import logging
logger = logging.getLogger()
logger.setLevel(logging.CRITICAL)

def main():

    parser = argparse.ArgumentParser(description="MEZE: MetalloEnZymE FF-builder for alchemistry")

    parser.add_argument("-if",
                        "--input-pdb-file",
                        dest="protein",
                        type=functions.file_exists,
                        required=True,
                        help="input pdb file for the metalloenzyme/protein")
    
    parser.add_argument("-g",
                        "--group-name",
                        dest="group_name",
                        help="group name for system, e.g. vim2/kpc2/ndm1/etc.; default taken from input filename",
                        default=None)

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
                        type=functions.path_exists,
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
                        type=float,
                        default=20) 
    
    parser.add_argument("-bs",
                        "--box-shape",
                        dest="box_shape",
                        help="box shape",
                        default="cubic") 

    parser.add_argument("-st",
                        "--sampling-time",
                        dest="sampling_time",
                         help="sampling time in nanoseconds",
                        type=float,
                        default=4) 
    
    parser.add_argument("-r",
                        "--repeats",
                        dest="repeats",
                        help="number of AFE repeat runs",
                        type=int, 
                        default=3)

    parser.add_argument("-min",
                        "--minimisation-steps",
                        dest="min_steps",
                        help="number of minimisation steps for equilibration stage",
                        type=int,
                        default=5000)
    
    parser.add_argument("-snvt",
                        "--short-nvt-runtime",
                        dest="short_nvt",
                        help="runtime in ps for short NVT equilibration",
                        type=float,
                        default=5) 
    
    parser.add_argument("-nvt",
                        "--nvt-runtime",
                        dest="nvt",
                        help="runtime in ps for NVT equilibration",
                        type=float,
                        default=50)

    parser.add_argument("-npt",
                        "--npt-runtime",
                        dest="npt",
                        help="runtime in ps for NPT equilibration",
                        type=float,
                        default=200)
    
    parser.add_argument("--em-step",
                        dest="emstep",
                        help="Step size for energy minimisation",
                        type=float,
                        default=0.001)
    
    parser.add_argument("--em-tolerance",
                        dest="emtol",
                        help="kJ mol-1 nm-1, Maximum force tolerance for energy minimisation",
                        type=float,
                        default=1000)

    arguments = parser.parse_args()

    network = Network.Network(workdir=arguments.working_directory,
                              ligand_path=arguments.ligand_directory,
                              group_name=arguments.group_name,
                              protein_file=arguments.protein,
                              protein_path=arguments.protein_directory,
                              water_model=arguments.water_model,
                              protein_ff=arguments.forcefield,
                              ligand_ff=arguments.ligand_forcefield,
                              ligand_charge=arguments.ligand_charge,
                              engine=arguments.engine,
                              sampling_time=arguments.sampling_time,
                              box_edges=arguments.box_edges,
                              box_shape=arguments.box_shape,
                              min_steps=arguments.min_steps,
                              short_nvt=arguments.short_nvt,
                              nvt=arguments.nvt,
                              npt=arguments.npt,
                              min_dt=arguments.emstep,
                              min_tol=arguments.emtol,
                              repeats=arguments.repeats)
    
    prepared_network = network.prepare_network()

    functions.write_slurm_script(template_file="02_add_water.sh", 
                                 path=prepared_network.afe_input_directory, 
                                 log_dir=prepared_network.log_directory,
                                 protocol_file=prepared_network.protocol_file)
    
    functions.write_slurm_script(template_file="03_heat_meze.sh", 
                                 path=prepared_network.afe_input_directory, 
                                 log_dir=prepared_network.log_directory,
                                 protocol_file=prepared_network.protocol_file)
    
    functions.write_slurm_script(template_file="04_meze.sh",
                                 path=prepared_network.afe_input_directory, 
                                 log_dir=prepared_network.log_directory,
                                 protocol_file=prepared_network.protocol_file)
     
    prepared_network.write_afe_run_script()

    functions.write_slurm_script(template_file="slurm_all.sh",
                                 path=prepared_network.afe_input_directory,
                                 log_dir=prepared_network.log_directory,
                                 protocol_file=prepared_network.protocol_file,
                                 extra_lines={"NUMBER_OF_LIGANDS": prepared_network.n_ligands,
                                              "ENGINE": prepared_network.md_engine,
                                              "AFE_PATH": prepared_network.afe_input_directory,
                                              "LIGANDS_DATA_FILE": prepared_network.ligands_dat_file,
                                              "TRANSFORMATIONS_DATA_FILE": prepared_network.network_file})
    
    print("\n")
    print("Successfully prepared meze network for AFE calculations.")
    print(f"Run scripts saved in: {prepared_network.afe_input_directory}")

if __name__ == "__main__":
    main()