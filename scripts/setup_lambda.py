import warnings
warnings.filterwarnings("ignore")
# import BioSimSpace as bss
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
n_transformations = len(transformation_folders)

for i in range(n_transformations):
    transformation_folder = transformation_folders[i]
    unbound_directory = transformation_folder + "/unbound/"
    bound_directory = transformation_folder + "/bound/"

    unbound_windows = sorted(glob.glob(unbound_directory+"*"))
    bound_windows = sorted(glob.glob(bound_directory+"*"))


unbound_paths = [transformation_folder + "/unbound/" for transformation_folder in transformation_folders]
bound_paths = [transformation_folder + "/bound/" for transformation_folder in transformation_folders]
unbound_transformation_folders = [sorted(glob.glob(unbound_path+"*")) for unbound_path in unbound_paths]
bound_transformation_folders = [sorted(glob.glob(bound_path+"*")) for bound_path in bound_paths]

minimisation_file = full_path + system_name + "/equilibration/bound/system_1/min/min.mdp"
nvt_file = full_path + system_name + "/equilibration/bound/system_1/nvt/nvt.mdp"
npt_file = full_path + system_name + "/equilibration/bound/system_1/npt/npt.mdp"

with open(minimisation_file) as min:
    minimisation_lines = min.readlines()
with open(nvt_file) as nvt:
    nvt_lines = nvt.readlines()
with open(npt_file) as npt:
    npt_lines = npt.readlines()

get_mdp_options = lambda index, lines: [line.split("=")[index].lstrip().rstrip().strip("\n") for line in lines]
min_keys, min_values = get_mdp_options(0, minimisation_lines), get_mdp_options(-1, minimisation_lines)
nvt_keys, nvt_values = get_mdp_options(0, nvt_lines), get_mdp_options(-1, nvt_lines)
npt_keys, npt_values = get_mdp_options(0, npt_lines), get_mdp_options(-1, npt_lines)

for folder in unbound_transformation_folders:
    for window in folder:
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
        
        # get mdp options for minimisation, nvt and mdp from equilibration folders

        # create dictionaries of each
        # write new mdp files with fn.lambda_mdps()

