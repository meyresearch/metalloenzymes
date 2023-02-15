import warnings
warnings.filterwarnings("ignore")
import BioSimSpace as bss 
import pandas as pd
import MDAnalysis as mda
import argparse
import os
import glob
import functions as fn


parser = argparse.ArgumentParser(description="solvate AFE systems")
parser.add_argument("system",
                     type=str,
                     help="system name; this is used to find the folder containing input files")

parser.add_argument("protein",
                    type=str,
                    help="filename for .prm7 and .rst7 files in ../system/inputs/protein/")


parser.add_argument("-c",
                    "--charge",
                    help="ligand net charge",
                    type=fn.check_charge,
                    default=0)

parser.add_argument("-fi",
                    "--fileinput",
                    help="filetype for ligand files",
                    type=str,
                    choices=["sdf", "mol2", "pdb"],
                    default="sdf")


parser.add_argument("-fo",
                    "--fileoutput",
                    help="output filetype",
                    type=str,
                    choices=["GROMACS", "AMBER"],
                    default="GROMACS")

arguments = parser.parse_args()
system_name = arguments.system
protein_name = arguments.protein
output_filetype = arguments.fileoutput
ligand_charge = arguments.charge

if output_filetype.upper() == "GROMACS":
    output_filetypes = ["Gro87", "GroTop"]
else:
    output_filetypes = ["PRM7", "RST7"]

full_path = os.getcwd() + "/"
if "scripts" in full_path:
    full_path = full_path.replace("/scripts/", "/")

protein_path = full_path + system_name + "/inputs/protein/"
protein_file = protein_path + protein_name
filetype = arguments.fileinput

afe_folder_path = full_path + system_name + "/afe/"
ligands_file = afe_folder_path + "ligands.dat"
protocol_file = afe_folder_path + "protocol.dat"
network_file = afe_folder_path + "network.dat"
ligand_path = full_path + system_name + "/inputs/ligands/"

# check files exist

with open(ligands_file) as file:
    ligand_lines = file.readlines()
ligand_names = [line.rstrip() for line in ligand_lines]
n_ligands = len(ligand_names)
with open(protocol_file) as file:
    protocol_lines = file.readlines()

ligand_force_field = protocol_lines[0].split()[-1]
solvent_force_field = protocol_lines[2].split()[-1]

box_edges = protocol_lines[3].split()[-1].split("*")[0]
box_axis_unit = bss.Units.Length.angstrom
box_type = protocol_lines[4].split()[-1]

ligand_files = glob.glob(f"{ligand_path}*.{filetype}")

protein = bss.IO.readMolecules([protein_file + ".rst7", protein_file + ".prm7"])[0]
ligands = [bss.IO.readMolecules(ligand_file)[0] for ligand_file in ligand_files]

for i in range(n_ligands):

    ligand_number = ligand_names[i].split("_")[1]
    print(f"working on ligand {ligand_number}")
    ligand_parameters = bss.Parameters.gaff2(ligands[i], net_charge=ligand_charge).getMolecule()

    box_min, box_max = ligand_parameters.getAxisAlignedBoundingBox()
    box_size = [y - x for x, y in zip(box_min, box_max)]
    box_sizes = [x + int(box_edges) * box_axis_unit for x in box_size] 

    system = ligand_parameters + protein

    system_box_min, system_box_max = system.getAxisAlignedBoundingBox()
    system_box_size = [y - x for x, y in zip(system_box_min, system_box_max)]
    system_box_sizes = [x + int(box_edges) * box_axis_unit for x in system_box_size]
    print("solvating unbound ligand")
    box, angles = bss.Box.cubic(max(box_sizes))
    ligand_solvated = bss.Solvent.solvate(solvent_force_field, molecule=ligand_parameters, box=box, angles=angles)
    box, angles = bss.Box.cubic(max(system_box_sizes))
    print("solvating bound system")
    system_solvated = bss.Solvent.solvate(solvent_force_field, molecule=system, box=box, angles=angles)

    ligand_savename = ligand_path + "ligand_" + ligand_number + "_solvated"
    system_savename = protein_path+ "system_" + ligand_number + "_solvated"
    bss.IO.saveMolecules(ligand_savename, ligand_solvated, output_filetypes)
    bss.IO.saveMolecules(system_savename, system_solvated, output_filetypes)    
