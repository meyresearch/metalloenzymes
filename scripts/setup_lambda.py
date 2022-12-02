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
                    type=fn.check_positive,
                    default=5000)

parser.add_argument("-nvt",
                    "--nvt-steps",
                    type=fn.check_positive,
                    help="number of steps for NVT equilibration with dt = 0.002 ps",
                    default=50000)

parser.add_argument("-npt",
                    "--npt-steps",
                    type=fn.check_positive,
                    help="number of steps for NPT equilibration with dt = 0.002 ps",
                    default=50000)

arguments = parser.parse_args()
system_name = arguments.system
minimisation_steps = arguments.minimisation_steps
nvt_steps = arguments.nvt_steps
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
    nvt_file = get_state_mdp(state, "nvt")
    npt_file = get_state_mdp(state, "npt")
    window_folders_list = [sorted(glob.glob(state_path+"*")) for state_path in state_paths]
    for path in window_folders_list:
        for window in path:
            minimisation_folder = window + "/min/"
            nvt_folder = window + "/nvt/"
            npt_folder = window + "/npt/"
            afe_folder = window + "/afe/"
            fn.create_dirs(minimisation_folder)
            fn.create_dirs(nvt_folder)
            fn.create_dirs(npt_folder)
            fn.create_dirs(afe_folder)    
            gromacs_files = glob.glob(window+"/gromacs*")
            for file in gromacs_files:
                try:
                    filename = file.split("/")[-1]
                    os.replace(file, afe_folder+filename)
                except FileNotFoundError:
                    pass
            afe_mdp_file = afe_folder+"gromacs.mdp"
            minimisation_options = fn.dict_from_mdp(minimisation_file)
            minimisation_options["nsteps"] = minimisation_steps
            nvt_options = fn.dict_from_mdp(nvt_file)
            nvt_options["nsteps"] = nvt_steps
            npt_options = fn.dict_from_mdp(npt_file)
            npt_options["nsteps"] = npt_steps
            fn.lambda_mdps(afe_mdp_file, minimisation_folder, "min", minimisation_options)   
            fn.lambda_mdps(afe_mdp_file, nvt_folder, "nvt", nvt_options)   
            fn.lambda_mdps(afe_mdp_file, npt_folder, "npt", npt_options)  
