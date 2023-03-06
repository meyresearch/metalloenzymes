import glob
import MDAnalysis as mda
import MDAnalysis.analysis.distances

files = glob.glob("../ligand_*/vim2_complex.amber.pdb")

zinc_oxygen_distances = []
for file in files:
    print(file.split("/")[-2])
    universe = mda.Universe(file)
    ligand = universe.select_atoms("resname LIG")
    zinc_1 = universe.select_atoms("resid 231").positions[0]
    zinc_2 = universe.select_atoms("resid 232").positions[0]

    o1 = ligand.select_atoms("id 3429").positions[0]
    o2 = ligand.select_atoms("id 3430").positions[0]
    o3 = ligand.select_atoms("id 3431").positions[0]

    get_distance = lambda zinc, oxygen: mda.analysis.distances.distance_array(zinc, oxygen)[0][0]
    zinc_1_distances = []
    zinc_1_distances.append(get_distance(zinc_1, o1))
    zinc_1_distances.append(get_distance(zinc_1, o2))
    zinc_1_distances.append(get_distance(zinc_1, o3))
    zinc_oxygen_distances.append(min(zinc_1_distances))

    zinc_2_distances = []
    zinc_2_distances.append(get_distance(zinc_2, o1))
    zinc_2_distances.append(get_distance(zinc_2, o2))
    zinc_2_distances.append(get_distance(zinc_2, o3))
    zinc_oxygen_distances.append(min(zinc_2_distances))
    print(zinc_oxygen_distances)
largest_zinc_oxygen_distance = max(zinc_oxygen_distances)
print(f"Largest zinc oxygen distance is: {largest_zinc_oxygen_distance:.3f} Å")
