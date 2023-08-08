#!/bin/bash


parent_dir=$1
protein_dir=$2

for i in {1..16}
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
	cat protein_only.pdb ZN1.pdb ZN2.pdb LIG.pdb WAT.pdb > vim2_cat.pdb
	pdb4amber vim2_cat.pdb > vim2_fixed.pdb
		
done

