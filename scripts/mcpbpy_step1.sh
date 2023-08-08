#!/bin/bash

export AMBERHOME="/home/jguven/Software/amber22/"

parent_dir=$1
protein_dir=$2

for i in {1..16}
do
    	ligand_dir="$parent_dir/ligand_$i"
	cd $ligand_dir
	cp $protein_dir/vim2.in .
	MCPB.py -i vim2.in -s 1
done
