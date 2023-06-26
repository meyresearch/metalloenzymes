import MDAnalysis as mda
import glob

paths = sorted(glob.glob("/home/jguven/projects/metalloenzymes/nonbonded_model_vim2/ligand_*/"))
# paths = ["/home/jguven/projects/metalloenzymes/nonbonded_model_vim2/ligand_8/"]
for path in paths:

    universe = mda.Universe(path + "vim2_solv.pdb")
    xtal_water = universe.select_atoms("resname HOH")
    water_oxygen_index = str(xtal_water.atoms.ids[0])

    ligand_oxygen = universe.select_atoms("resname L8J and name O13")
    ligand_oxygen_index = str(ligand_oxygen.atoms.ids[0])

    with open(path + "disres.rst", "r") as file:
        lines = file.readlines()
    
    water_coordination = lines[3].split(",")[1]
    ligand_coordination = lines[-1].split(",")[1]
    lines[3] = lines[3].replace(water_coordination, water_oxygen_index)
    lines[-1] = lines[-1].replace(ligand_coordination, ligand_oxygen_index)

    with open(path + "disres.rst", "w") as file:
        file.writelines(lines)