import BioSimSpace as bss


protein = bss.IO.readPDB("../protein/vim2.amber.pdb")[0]
zn1 = bss.IO.readPDB("../protein/ZN1.pdb")[0]
zn2 = bss.IO.readPDB("../protein/ZN2.pdb")[0]
water = bss.IO.readPDB("../protein/WAT.pdb")[0]

complex = protein + zn1 + zn2 + water

bss.IO.saveMolecules("../protein/complex.pdb", complex, "pdb")
