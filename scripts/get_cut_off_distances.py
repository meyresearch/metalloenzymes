import glob
import MDAnalysis as mda
import MDAnalysis.analysis.distances


files = glob.glob("../dft/ligand_*/vim2_complex.amber.pdb")
for file in files:
    print(file)
    universe = mda.Universe(file)
    ligand = universe.select_atoms("resname LIG")
    zinc_1 = universe.select_atoms("resid 231").positions[0]
    zinc_2 = universe.select_atoms("resid 232").positions[0]
    
    o1 = ligand.select_atoms("id 3429").positions[0]
    o2 = ligand.select_atoms("id 3430").positions[0]
    o3 = ligand.select_atoms("id 3431").positions[0]

    get_distance = lambda zinc, oxygen: mda.analysis.distances.distance_array(zinc, oxygen)

    print("zinc 1")
    print(get_distance(zinc_1, o1))
    print(get_distance(zinc_1, o2))
    print(get_distance(zinc_1, o3))
    print("zinc 2")
    print(get_distance(zinc_2, o1))
    print(get_distance(zinc_2, o2))
    print(get_distance(zinc_2, o3))

    
