from qmMeze import Meze
import argparse
import functions
import os




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
    
    parser.add_argument("-wf",
                        "--water-file",
                        dest="water_file",
                        help="water file, default is water.pdb",
                        default=None)
    
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


    parser.add_argument("-min",
                        "--minimisation-steps",
                        dest="min_steps",
                        help="number of minimisation steps for equilibration stage",
                        type=int,
                        default=1000)
    
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
    
    parser.add_argument("-m",
                        "--metal",
                        dest="metal", 
                        help="name of the metal ion",
                        default="ZN",
                        type=str)
    
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

    arguments = parser.parse_args()


    meze = Meze(metal=arguments.metal,
                is_qm=True,
                cut_off=arguments.cut_off,
                force_constant_0=arguments.force_constant_0,
                workdir=arguments.working_directory,
                ligand_path=arguments.ligand_directory,
                group_name=arguments.group_name,
                protein_file=arguments.protein,
                protein_path=arguments.protein_directory,
                water_model=arguments.water_model,
                protein_ff=arguments.forcefield,
                ligand_ff=arguments.ligand_forcefield,
                ligand_charge=arguments.ligand_charge,
                sampling_time=arguments.sampling_time,
                box_edges=arguments.box_edges,
                box_shape=arguments.box_shape,
                min_steps=arguments.min_steps,
                short_nvt=arguments.short_nvt,
                nvt=arguments.nvt,
                npt=arguments.npt,
                min_dt=arguments.emstep,
                min_tol=arguments.emtol,
                water_file=arguments.water_file)

    prepared_meze = meze.prepare_meze()
    
    functions.write_slurm_script(template_file="02_add_water.sh", 
                                 path=prepared_meze.input_directory, 
                                 log_dir=prepared_meze.log_directory,
                                 protocol_file=prepared_meze.protocol_file)

    
    ### Write slurm_all.sh script



if __name__ == "__main__":
   main()