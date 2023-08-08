#!/bin/bash

ligands=(  1 2 3 4 5 6 7 8  9 10 11 12 13 14 15 16 )

parent_dir="$HOME/projects/qmmm/nonbonded_models/integer_charge/part_protonated/"
protein_dir="$HOME/projects/qmmm/nonbonded_models/integer_charge/part_protonated/protein"

for i in {1..16}
do
	ligand_directory="$parent_dir/ligand_$i"
	cd $ligand_directory

	cp $protein_dir/protein_only.pdb .
	cp $protein_dir/WT1.pdb .
	cp $protein_dir/ZN.pdb .

	antechamber -fi pdb -fo mol2 -i WT1.pdb -o WT1.mol2 -at amber -c bcc -pf y 
	sed -i "s/HO/HW/g" WT1.mol2
	metalpdb2mol2.py -i ZN.pdb -o ZN.mol2 -c 2
	# mv "ligand_${i}.pdb" L8J.pdb
	cat protein_only.pdb ZN.pdb L8J.pdb WT1.pdb > vim2_cat.pdb
	pdb4amber vim2_cat.pdb > vim2_fixed.pdb
		
done

