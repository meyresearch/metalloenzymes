import warnings
warnings.filterwarnings("ignore")
import BioSimSpace as bss
import argparse
import os
import sys
import glob
import csv
import pandas as pd


def create_complex(protein_path: str, protein_file: str, water_file: str, output=None) -> str:
    protein = bss.IO.readMolecules(protein_path+protein_file)
    if water_file is not None and output == "":
        output = "protein_complex"
        xtal_water = bss.IO.readMolecules(protein_path+water_file)
        protein_complex = protein + xtal_water
        bss.IO.saveMolecules(protein_path+output, protein_complex, fileformat="pdb")
    elif water_file is not None and output != "":
        xtal_water = bss.IO.readMolecules(protein_path+water_file)
        protein_complex = protein + xtal_water
        bss.IO.saveMolecules(protein_path+output, protein_complex, fileformat="pdb")
    elif water_file is None and output == "":
        output = "protein_complex"
        protein_complex = protein
        bss.IO.saveMolecules(output, protein_complex, fileformat="pdb")
    elif water_file is None and output != "":
        protein_complex = protein
        bss.IO.saveMolecules(protein_path+output, protein_complex, fileformat="pdb")
    return output


def is_file(file: str) -> None:
    if not os.path.isfile(file):
        print(f"The file {file} does not exist")
        sys.exit()


parser = argparse.ArgumentParser(description="prepare AFE calculations")

parser.add_argument("system",
                     type=str,
                     help="system name; this is used to find the folder containing input files")

parser.add_argument("protein",
                    type=str,
                    help="protein pdb file in ../system/inputs/protein/")

parser.add_argument("-w",
                    "--water",
                    type=str,
                    help="water pdb file in ../system/inputs/protein/")

parser.add_argument("-ff",
                    "--forcefield",
                    type=str,
                    help="forcefield for protein",
                    choices=["ff14SB", "zaff"],
                    default="ff14SB")

parser.add_argument("-s",
                    "--solvent",
                    type=str,
                    help="solvent",
                    default="tip3p")

parser.add_argument("-o",
                    "--output",
                    type=str,
                    help="output file",
                    default="")

arguments = parser.parse_args()
system = arguments.system
protein_path = "../" + system + "/inputs/protein/"
ligand_path = system + "/inputs/ligands/"
protein_file = arguments.protein
water_file = arguments.water
output_file = arguments.output
forcefield = arguments.forcefield
solvent = arguments.solvent.lower()
tleap_file = protein_path + "tleap.in"
is_file(protein_path+protein_file)
complex_file = create_complex(protein_path, protein_file, water_file, output_file)

if arguments.water is not None:
    water_file = protein_path + arguments.water
    is_file(water_file)
else: 
    water_file = None

if forcefield.lower() == "zaff":
    try:
        tleap_command = f"tleap -s -f {tleap_file} > {protein_path}" + "tleap.out"
        os.system(tleap_command)
    except FileNotFoundError:
        print("tleap input file for zaff does not exist and is required with --forcefield=zaff")
        sys.exit()
elif forcefield.lower() == "ff14sb":
    forcefield = "ff14SB"
    with open(tleap_file, "w") as tleap_in:
        tleap_in.write(f"source leaprc.protein.{forcefield}\n")
        tleap_in.write(f"source leaprc.water.{solvent}\n")
        tleap_in.write(f"complex = loadpdb {protein_path}{complex_file}.pdb\n")
        tleap_in.write(f"saveamberparm complex {protein_path+system}_tleap.prm7 {protein_path+system}_tleap.rst7\n")
        tleap_in.write("quit")
    tleap_command = f"tleap -s -f {tleap_file} > {protein_path}" + "tleap.out"
    os.system(tleap_command)

os.chdir("../")
ligand_files = sorted(glob.glob(f"{ligand_path}docked_*.sdf"))
ligands = [bss.IO.readMolecules(filepath)[0] for filepath in ligand_files]
ligand_names = [filepath.split("/")[-1].replace(".sdf","") for filepath in ligand_files]
transformations, lomap_scores = bss.Align.generateNetwork(ligands, plot_network=True, names=ligand_names, work_dir=ligand_path)

pert_network_dict = {}
transformations_named = [(ligand_names[transf[0]], ligand_names[transf[1]]) for transf in transformations]
for transf, score in zip(transformations_named, lomap_scores):
    pert_network_dict[transf] = score

with open(ligand_path+f"/lomap_{system}.csv", "w") as lomap_out:
    for key, value in pert_network_dict.items():
        lomap_out.write(f"{key}: {value}\n")
