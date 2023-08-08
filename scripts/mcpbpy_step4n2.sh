#!/bin/bash

export AMBERHOME="/home/jguven/Software/amber22/"

ligands=( 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 )
parent_dir="$HOME/projects/qmmm/nonbonded_models/integer_charge/part_protonated/"
protein_dir="$HOME/projects/qmmm/nonbonded_models/integer_charge/part_protonated/protein/"

for i in {1..16}
do
    	ligand_dir="$parent_dir/ligand_$i"
    	cd $ligand_dir
	cp $protein_dir/vim2.in .
	MCPB.py -i vim2.in -s 4n2
done
