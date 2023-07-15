#!/bin/bash

ligands=(  1 2 3 4 5 6 7 8  9 10 11 12 13 14 15 16 )

parent_dir="$HOME/projects/alchemistry/vim2_deprotonated_ligand/parameterisation"
protein_dir="$HOME/projects/alchemistry/vim2_deprotonated_ligand/inputs/protein"
# prep_dir="$parent_dir/system_preparation"
# steps=( 1 2 )

for i in ${ligands[@]}
do
	ligand_directory="$parent_dir/ligand_$i"
	cd $ligand_directory

	cp $protein_dir/protein_only.pdb .
	cp $protein_dir/WAT.pdb .
	cp $protein_dir/ZN*.pdb .

	antechamber -fi pdb -fo mol2 -i WAT.pdb -o WAT.mol2 -at amber -c bcc -pf y 
	sed -i "s/HO/HW/g" WAT.mol2
	metalpdb2mol2.py -i ZN1.pdb -o ZN1.mol2 -c 2
	metalpdb2mol2.py -i ZN2.pdb -o ZN2.mol2 -c 2
	# mv "ligand_${i}.pdb" L8J.pdb
	cat protein_only.pdb ZN1.pdb ZN2.pdb L8J.pdb WAT.pdb > vim2_cat.pdb
	pdb4amber vim2_cat.pdb > vim2_fixed.pdb
		
done

