import sys
import os


lig = sys.argv[1]
full_path = os.getenv("HOME") + "/projects/metalloenzymes/dft"

with open(f"{full_path}/ligand_{lig}/01_minimisation/vim2_solv.top", "r") as topfile:
    lines = topfile.readlines()

info = [line.split() for line in lines]
sections = []
for i in range(len(info)):
    if len(info[i]) > 0:
        if "moleculetype" in info[i]:
            sections.append(i)

position_restraints = ["; Include Position restraint file\n",
                       "#ifdef POSRES\n",
                       "#include \"../02_nvt/posre.itp\"\n",
                       "#endif\n",
                       "\n"]


with open(f"{full_path}/ligand_{lig}/01_minimisation/vim2_solv.top", "w") as topfile:
    for row in lines[:sections[1]]:
        topfile.writelines(row)
    for line in position_restraints:
        topfile.writelines(line)
    for row in lines[sections[1]:]:
        topfile.writelines(row)
