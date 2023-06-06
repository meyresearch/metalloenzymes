import argparse
import functions as fn
import os
import glob


parser = argparse.ArgumentParser(description="setup minimisation and equilibration in each lambda window")

parser.add_argument("system",
                     type=str,
                     help="system name; this is used to find the folder containing input files")

parser.add_argument("-m",
                    "--minimisation-steps",
                    type=fn.check_positive_integer,
                    default=5000)

parser.add_argument("-e",
                    "--em-step",
                    help="stepsize for energy minimisation",
                    type=fn.check_positive_float,
                    default=0.01)

parser.add_argument("-f",
                    "--f-max",
                    help="maximum force, emtol in mdp options",
                    type=fn.check_positive_float,
                    default=1000)

parser.add_argument("-nvt",
                    "--nvt-steps",
                    type=fn.check_positive_integer,
                    help="number of steps for NVT equilibration with dt = 0.002 ps",
                    default=50000)

parser.add_argument("-t",
                    "--time-step",
                    type=fn.check_positive_float,
                    help="NVT timestep for NVT equilibration (dt)",
                    default=0.0002)

parser.add_argument("-npt",
                    "--npt-steps",
                    type=fn.check_positive_integer,
                    help="number of steps for NPT equilibration with dt = 0.002 ps",
                    default=50000)

arguments = parser.parse_args()
system_name = arguments.system
minimisation_steps = arguments.minimisation_steps
nvt_steps = arguments.nvt_steps
dt = arguments.time_step
em_step = arguments.em_step
em_tolerance = arguments.f_max
npt_steps = arguments.npt_steps

full_path = os.getcwd() + "/"
if "scripts" in full_path:
    full_path = full_path.replace("/scripts/", "/")

afe_folder_path = full_path + system_name + "/afe/"
protocol_file = afe_folder_path + "protocol.dat"
with open(protocol_file) as protocol:
    protocol_lines = protocol.readlines()

settings = [line.split("=")[-1].strip() for line in protocol_lines]
engine = settings[-1]
gmx_output_path = full_path + system_name + "/outputs/" + engine.upper() + "/"

transformation_folders = glob.glob(gmx_output_path+"*")
states = ["unbound", "bound"]
get_state_mdp = lambda state, process: full_path + system_name + f"/equilibration/{state}/{process}.mdp"

for state in states:
    state_paths = [transformation_folder + "/" + state + "/" for transformation_folder in transformation_folders]
    minimisation_file = get_state_mdp(state, "min")
    r_nvt_file = get_state_mdp(state, "r_nvt")
    nvt_file = get_state_mdp(state, "nvt")
    npt_file = get_state_mdp(state, "npt")
    window_folders_list = [sorted(glob.glob(state_path+"*")) for state_path in state_paths]
    for path in window_folders_list:
        for window in path:
            minimisation_folder = window + "/min/"
            r_nvt_folder = window + "/r_nvt/"
            nvt_folder = window + "/nvt/"
            npt_folder = window + "/npt/"
            afe_folder = window + "/afe/"
            fn.create_dirs(minimisation_folder)
            fn.create_dirs(r_nvt_folder)
            fn.create_dirs(nvt_folder)
            fn.create_dirs(npt_folder)
            fn.create_dirs(afe_folder)    
            gromacs_files = glob.glob(window+"/gromacs*")
            try:
                os.replace(window+"/gromacs.gro", minimisation_folder+"/min.gro")
                os.replace(window+"/gromacs.top", minimisation_folder+"/min.top")
            except FileNotFoundError as error:
                pass
            for file in gromacs_files:
                try:
                    filename = file.split("/")[-1]
                    os.replace(file, afe_folder+filename)
                except FileNotFoundError:
                    pass
            afe_mdp_file = afe_folder+"gromacs.mdp"
            minimisation_options = fn.dict_from_mdp(minimisation_file)
            minimisation_options["nsteps"] = minimisation_steps
            minimisation_options["emstep"] = em_step 
            minimisation_options["emtol"] = em_tolerance  
            r_nvt_options = fn.dict_from_mdp(r_nvt_file)
            r_nvt_options["dt"] = dt
            nvt_options = fn.dict_from_mdp(nvt_file)
            nvt_options["nsteps"] = nvt_steps
            npt_options = fn.dict_from_mdp(npt_file)
            npt_options["nsteps"] = npt_steps
            fn.lambda_mdps(afe_mdp_file, minimisation_folder, "min", minimisation_options)   
            fn.lambda_mdps(afe_mdp_file, r_nvt_folder, "r_nvt", r_nvt_options)    
            fn.lambda_mdps(afe_mdp_file, nvt_folder, "nvt", nvt_options)   
            fn.lambda_mdps(afe_mdp_file, npt_folder, "npt", npt_options)  

# bb restrained nvt for bound
bound_paths = [transformation_folder + "/" + "bound" + "/" for transformation_folder in transformation_folders]
# min2_file = get_state_mdp("bound", "min")
bb_r_nvt_file = get_state_mdp("bound", "bb_r_nvt")
window_folders_list = [sorted(glob.glob(bound_path+"*")) for bound_path in bound_paths]
for path in window_folders_list:
    for window in path:
        afe_folder = window + "/afe/"
        afe_mdp_file = afe_folder+"gromacs.mdp"
        # min2_folder = window + "/min_2/"
        bb_r_nvt_folder = window + "/bb_r_nvt/"
        # fn.create_dirs(min2_folder)    
        fn.create_dirs(bb_r_nvt_folder)
        # min2_options = fn.dict_from_mdp(min2_file)
        # min2_options["nsteps"] = 20000
        # min2_options["emstep"] = em_step 
        bb_r_nvt_options = fn.dict_from_mdp(bb_r_nvt_file)
        fn.lambda_mdps(afe_mdp_file, bb_r_nvt_folder, "bb_r_nvt", bb_r_nvt_options)   
        # fn.lambda_mdps(afe_mdp_file, min2_folder, "min_2", min2_options)
