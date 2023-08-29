# MetalloEnZymE parameterisation program (MEZE)

**Authors**: J. Jasmin Güven


# 1. System setup

- Start with a pdb file of the protein. You can save it in e.g. `inputs/protein/` in your project directory
- If your protein contains metal ions, save these in pdb format in the above directory
- Create ligand files in `*.sdf` or `*.mol2` formats and save in the same directory, e.g. `inputs/ligands/`
- If you want to keep crystallographic waters, save these in the protein directory in pdb format

# 2. Prepare meze

- Run `meze.py` to prepare your protein for AFE calculations