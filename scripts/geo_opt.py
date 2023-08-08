import psi4

file = "/home/jguven/projects/parameterisation_with_psi4/vim2_small_opt.xyz"
with open(file, "r") as f:
    xyz_contents = f.read()

geometry = psi4.core.Molecule.from_string(xyz_contents, dtype="xyz")
charge = geometry.molecular_charge()
multiplicity = geometry.multiplicity()
print(f"Charge is {charge}")
print(f"Multiplicity is {multiplicity}") 

psi4.set_memory("12 GB")
psi4.set_num_threads(30)

psi4.optimize(name="B3LYP/6-31G*", molecule=geometry)