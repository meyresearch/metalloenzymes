import BioSimSpace as bss

parent_dir = "/home/jguven/projects/zn_models/nonbonded_models/deprotonated/charge_fitting/"
protein = bss.IO.readMolecules(f"{[parent_dir]}/vim2_hpp_output.pdb")[0]
zn1 = bss.IO.readMolecules(f"{parent_dir}/ZN1.pdb")[0]
zn2 = bss.IO.readMolecules(f"{parent_dir}/ZN2.pdb")[0]
water = bss.IO.readMolecules(f"{parent_dir}/WAT.pdb")[0]

complex = protein + zn1 + zn2 + water

bss.IO.saveMolecules("../protein/bss_complex.pdb", complex, "pdb")
