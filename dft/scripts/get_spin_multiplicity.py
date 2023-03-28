import mendeleev 
import MDAnalysis as mda
import argparse 


parser = argparse.ArgumentParser("calculate the number of electrons in a pdb file")
parser.add_argument("-fi", "--file-input", type=str, help="full path to pdb file")
parser.add_argument("-zn", "--zinc-oxidation", type=int, help="oxidation state for zinc ions", default=2)
parser.add_argument("-c", "--charge", type=int, help="total formal charge in the system")
arguments = parser.parse_args()



pdb_file = arguments.file_input
zinc_oxidation = arguments.zinc_oxidation
charge = arguments.charge

universe = mda.Universe(pdb_file)
elements = list(universe.atoms.types)

counts = {element: elements.count(element) for element in elements}


list_of_atomic_numbers = []
for element in counts:
    if len(element) == 2:
        second_letter = element[-1].lower()
        element = element[0] + second_letter
    list_of_atomic_numbers.append(mendeleev.element(element).atomic_number)

atomic_numbers = {list(counts.keys())[i]:list_of_atomic_numbers[i] for i in range(len(counts))}
electrons_per_element = {element: counts[element] * atomic_numbers[element] for element in counts}

total_n_electrons = sum(electrons_per_element.values())
number_of_zincs = counts["ZN"]

if charge < 0:
    final_n_electrons = total_n_electrons - zinc_oxidation*number_of_zincs + charge
elif charge > 0:
    final_n_electrons = total_n_electrons - zinc_oxidation*number_of_zincs - charge # Double-check
else:
    final_n_electrons = total_n_electrons - zinc_oxidation*number_of_zincs

print(f"Electrons: {total_n_electrons}, Zincs: {number_of_zincs}. Total n:o electrons: {final_n_electrons}")
