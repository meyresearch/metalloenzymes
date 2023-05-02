import BioSimSpace as bss


protein = bss.IO.readMolecules("../protein/vim2_amber_output.pdb")[0]
zn1 = bss.IO.readMolecules("../protein/ZN1.pdb")[0]
zn2 = bss.IO.readMolecules("../protein/ZN2.pdb")[0]
water = bss.IO.readMolecules("../protein/WAT.pdb")[0]

complex = protein + zn1 + zn2 + water

bss.IO.saveMolecules("../protein/bss_complex.pdb", complex, "pdb")
