import argparse
import functions
import os
import sofra
import meze
import logging
logger = logging.getLogger()
logger.setLevel(logging.CRITICAL)

def main():

    parser = argparse.ArgumentParser(description="MEZE: MetalloEnZymE FF-builder for alchemistry")

    parser.add_argument("-n",
                        "--non-metal",
                        dest="non_metal", 
                        help="boolean flag for preparing non-metal proteins",
                        action="store_true")
    
    parser.add_argument("-co",
                        "--cut-off",
                        dest="cut_off",
                        help="cut off distance in Å, used to define zinc coordinating ligands and active site",
                        default=2.6,
                        type=float)
    
    parser.add_argument("-fc0",
                        "--force-constant-0",
                        dest="force_constant_0",
                        help="force constant for model 0",
                        default=100,
                        type=float)    
    
    parser.add_argument("-if",
                        "--input-pdb-file",
                        dest="protein",
                        type=functions.file_exists,
                        help="input pdb file for the metalloenzyme/protein")
    
    parser.add_argument("-g",
                        "--group-name",
                        dest="group_name",
                        help="group name for system, e.g. vim2/kpc2/ndm1/etc.; default taken from input filename",
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
    
    parser.add_argument("-lw",
                        "--lambda-windows",
                        dest="lambdas",
                        help="number of lambda windows to use for AFE transformations",
                        type=int,
                        default=11)    
    
    parser.add_argument("-ld",
                        "--lambda-windows-difficult",
                        dest="n_difficult",
                        help="number of lambda windows to use for difficult AFE transformations",
                        type=int,
                        default=17) 
    
    parser.add_argument("-sm",
                        "--solvation-method",
                        dest="solvation_method",
                        help="MD engine used for solvation; default is BioSimSpace (GROMACS)",
                        choices=["amber", "gromacs"],
                        default="amber",
                        type=str)
    
    parser.add_argument("-sc",
                        "--solvent-closeness",
                        dest="solvent_closeness",
                        help="control how close, in \AA, solvent atoms can come to solute atoms, default 1.0",
                        type=float, 
                        default=1.0)

    parser.add_argument("-et",
                        "--extra-transformations",
                        dest="extra_transformations_file",
                        help="file containing additional transformations in format: lig A lig B; equivalent to BioSimSpace.generateNetwork links file",
                        type=str)
    
    parser.add_argument("-pf",
                        "--protocol-file",
                        dest="protocol_file",
                        help="protocol file",
                        type=str,
                        default=os.getcwd() + "/afe/protocol.dat")

    parser.add_argument("-cs",
                        "--cutoff-scheme",
                        dest="cutoff_scheme",
                        help="PME or Reaction Field",
                        default="rf",
                        choices=["pme", "rf"],
                        type=str)
    
    parser.add_argument("-os",
                        "--only-save-end-states",
                        help="tell somd to only save trajectories for lambda 0 and lambda 1",
                        action=argparse.BooleanOptionalAction,
                        dest="only_save_end_states")
    

    arguments = parser.parse_args()

    if not arguments.non_metal:
        metal = True
    elif arguments.non_metal:
        metal = False

    if not arguments.extra_transformations_file:
        if metal:
            network = meze.Meze(protein_file=arguments.protein,
                                cut_off=arguments.cut_off,
                                force_constant_0=arguments.force_constant_0,
                                workdir=arguments.working_directory,
                                ligand_path=arguments.ligand_directory,
                                ligand_charge=arguments.ligand_charge,
                                ligand_ff=arguments.ligand_forcefield,
                                group_name=arguments.group_name,        
                                protein_path=arguments.protein_directory,
                                water_model=arguments.water_model,
                                protein_ff=arguments.forcefield,
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
                                repeats=arguments.repeats,
                                n_normal=arguments.lambdas,
                                n_difficult=arguments.n_difficult,
                                cutoff_scheme=arguments.cutoff_scheme,
                                solvation_method=arguments.solvation_method,
                                solvent_closeness=arguments.solvent_closeness, 
                                only_save_end_states=arguments.only_save_end_states)

        
        elif not metal:

            network = sofra.Sofra(protein_file=arguments.protein,
                                    workdir=arguments.working_directory,
                                    ligand_path=arguments.ligand_directory,
                                    ligand_charge=arguments.ligand_charge,
                                    ligand_ff=arguments.ligand_forcefield,
                                    group_name=arguments.group_name,        
                                    protein_path=arguments.protein_directory,
                                    water_model=arguments.water_model,
                                    protein_ff=arguments.forcefield,
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
                                    repeats=arguments.repeats,
                                    n_normal=arguments.lambdas,
                                    n_difficult=arguments.n_difficult,
                                    solvation_method=arguments.solvation_method,
                                    solvent_closeness=arguments.solvent_closeness, 
                                    only_save_end_states=arguments.only_save_end_states)

        prepared_network = network.prepare_network()

        functions.write_slurm_script(template_file="02_add_water.sh", 
                                     path=prepared_network.afe_input_directory, 
                                     log_dir=prepared_network.log_directory,
                                     protocol_file=prepared_network.protocol_file,
                                     extra_lines={"LIGANDS_DATA_FILE": prepared_network.ligands_dat_file,
                                                  "TRANSFORMATIONS_DATA_FILE": prepared_network.network_file})
        
        functions.write_slurm_script(template_file="03_heat_meze.sh", 
                                    path=prepared_network.afe_input_directory, 
                                    log_dir=prepared_network.log_directory,
                                    protocol_file=prepared_network.protocol_file,
                                    extra_lines={"LIGANDS_DATA_FILE": prepared_network.ligands_dat_file,
                                                 "TRANSFORMATIONS_DATA_FILE": prepared_network.network_file})
        
        functions.write_slurm_script(template_file="04_meze.sh",
                                    path=prepared_network.afe_input_directory, 
                                    log_dir=prepared_network.log_directory,
                                    protocol_file=prepared_network.protocol_file)
        
        prepared_network.write_afe_run_script()

        functions.write_slurm_script(template_file="06_analyse.sh",
                                    path=prepared_network.afe_input_directory,
                                    log_dir=prepared_network.log_directory,
                                    protocol_file=prepared_network.protein_file)

        functions.write_slurm_script(template_file="slurm_all.sh",
                                    path=prepared_network.afe_input_directory,
                                    log_dir=prepared_network.log_directory,
                                    protocol_file=prepared_network.protocol_file,
                                    extra_lines={"NUMBER_OF_LIGANDS": prepared_network.n_ligands,
                                                 "ENGINE": prepared_network.md_engine,
                                                 "AFE_PATH": prepared_network.afe_input_directory,
                                                 "LIGANDS_DATA_FILE": prepared_network.ligands_dat_file,
                                                 "TRANSFORMATIONS_DATA_FILE": prepared_network.network_file})
    
    elif arguments.extra_transformations_file:

        protocol_file = functions.file_exists(arguments.protocol_file)
        protocol = functions.input_to_dict(protocol_file)
        if metal:
            prepared_network = meze.Meze(prepared=True,
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
                                         n_normal=arguments.lambdas,
                                         n_difficult=arguments.n_difficult, 
                                         only_save_end_states=arguments.only_save_end_states)
            
        elif not metal:
            prepared_network = sofra.Sofra(prepared=True,
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
                                           pressure=protocol["pressure"],
                                           n_normal=arguments.lambdas,
                                           n_difficult=arguments.n_difficult,
                                           solvation_method=protocol["solvation_method"],
                                           solvent_closeness=protocol["solvent_closeness"],
                                           n_difficult=arguments.n_difficult, 
                                           only_save_end_states=arguments.only_save_end_states)

        links_file = functions.file_exists(arguments.extra_transformations_file)

        _, transformations_file = prepared_network.set_transformations(links_file=links_file)

        functions.write_slurm_script(template_file="run_extra_edges.sh",
                                     path=prepared_network.afe_input_directory,
                                     log_dir=prepared_network.log_directory,
                                     protocol_file=protocol_file,
                                     extra_lines={"ENGINE": prepared_network.md_engine,
                                                  "AFE_PATH": prepared_network.afe_input_directory,
                                                  "TRANSFORMATIONS_DATA_FILE": transformations_file})  
          
        functions.write_slurm_script(template_file="07_extra_meze.sh",
                                     path=prepared_network.afe_input_directory, 
                                     log_dir=prepared_network.log_directory,
                                     protocol_file=arguments.protocol_file,
                                     extra_lines={"EXTRA_TRANSFORMATIONS": transformations_file})
    print("\n")
    print("Successfully prepared meze network for AFE calculations.")
    print(f"Run scripts saved in: {prepared_network.afe_input_directory}")

if __name__ == "__main__":
    main()