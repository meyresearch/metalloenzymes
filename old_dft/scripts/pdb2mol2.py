import BioSimSpace as bss
import argparse

parser = argparse.ArgumentParser(description="open renamed antechamber output file save as mol2")

parser.add_argument("pdb_file",
                    type=str,
                    help="input pdb file to be converted to mol2")

arguments = parser.parse_args()
input = arguments.pdb_file
system = bss.IO.readMolecules(input)
bss.IO.saveMolecules("LIG", system, "mol2")
