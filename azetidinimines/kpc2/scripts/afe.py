import argparse
import os
import glob
import pandas as pd



parser = argparse.ArgumentParser(description="setup minimisation and equilibration in each lambda window")

parser.add_argument("system",
                     type=str,
                     help="system name; this is used to find the folder containing input files")
parser.add_argument("network",
                    type=str,
                    help="filename including full path to network file")

arguments = parser.parse_args()
system_name = arguments.system
network_file = arguments.network


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

transformation_paths = glob.glob(gmx_output_path+"*")

unbound_path = [transformation_folder + "/unbound/" for transformation_folder in transformation_paths]
bound_path = [transformation_folder + "/bound/" for transformation_folder in transformation_paths]

unbound_windows = [sorted(glob.glob(state_path+"*")) for state_path in unbound_path]
bound_windows = [sorted(glob.glob(state_path+"*")) for state_path in bound_path]

afe_runs = afe_folder_path + "run_unbound.dat"
with open(afe_runs, "w") as file:
    for i in range(len(unbound_windows)):
        for unbound_line in unbound_windows[i]:
            file.write(unbound_line+"\n")


afe_runs = afe_folder_path + "run_bound.dat"
with open(afe_runs, "w") as file:
    for i in range(len(bound_windows)):
        for bound_line in bound_windows[i]:
            file.write(bound_line+"\n")
