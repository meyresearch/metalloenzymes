import warnings
warnings.filterwarnings("ignore")
import BioSimSpace as bss
import argparse
import os
import sys
import glob


def create_complex(protein_file: str, water_file: str, output="protein_complex") -> str:
    protein = bss.IO.readMolecules(protein_file)
    xtal_water = bss.IO.readMolecules(water_file)
    protein_complex = protein + xtal_water
    bss.IO.saveMolecules(output, protein_complex, fileformat="pdb")
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

parser.add_argument("water",
                    type=str,
                    help="water pdb file in ../system/inputs/protein/")

parser.add_argument("-ff",
                    "--forcefield",
                    type=str,
                    help="forcefield for protein",
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
# ligand_path = "../" + system + "/inputs/ligands/"
ligand_path = system + "/inputs/ligands/"
protein_file = protein_path + arguments.protein
water_file = protein_path + arguments.water
output_file = protein_path + arguments.output
forcefield = arguments.forcefield
solvent = arguments.solvent.lower()
tleap_file = protein_path + "tleap.in"
is_file(protein_file)
is_file(water_file)

if output_file != "":
    complex_file = create_complex(protein_file, water_file, output_file)
else:
    complex_file = create_complex(protein_file, water_file)

with open(tleap_file, "w") as tleap_in:
    tleap_in.write(f"source leaprc.protein.{forcefield}\n")
    tleap_in.write(f"source leaprc.water.{solvent}\n")
    tleap_in.write(f"complex = loadpdb {complex_file}.pdb\n")
    tleap_in.write(f"saveamberparm complex {protein_path+system}_tleap.prm7 {protein_path+system}_tleap.rst7\n")
    tleap_in.write("quit")
tleap_command = f"tleap -s -f {tleap_file} > {protein_path}" + "tleap.out"
os.system(tleap_command)
# print(os.getcwd())
os.chdir("../")
ligand_files = sorted(glob.glob(f"{ligand_path}docked_*.sdf"))
ligands = [bss.IO.readMolecules(filepath)[0] for filepath in ligand_files]
ligand_names = [filepath.split("/")[-1].replace(".sdf","") for filepath in ligand_files]
transformations, lomap_scores = bss.Align.generateNetwork(ligands, plot_network=True, names=ligand_names, work_dir=ligand_path)
pert_network_dict = {}
transformations_named = [(ligand_names[transf[0]], ligand_names[transf[1]]) for transf in transformations]
for transf, score in zip(transformations_named, lomap_scores):
    print(transf, score)
    pert_network_dict[transf] = score