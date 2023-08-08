#!/bin/bash

export AMBERHOME=/home/jguven/Software/amber22/

parent_dir=$HOME/projects/qmmm/nonbonded_models/integer_charge/part_protonated/


for i in {1..16}
do
	ligand_dir=$parent_dir/ligand_$i
	cd $ligand_dir
	cp vim2_tleap.in mcpbpy_vim2_tleap.in
	sed -i "3i source leaprc.gaff2" vim2_tleap.in
	sed -i "s/WT1.mol2/WT1/g" vim2_tleap.in
	sed -i "s/L8J.mol2/L8J/g" vim2_tleap.in
	tleap -s -f vim2_tleap.in > vim2_tleap.out
done

