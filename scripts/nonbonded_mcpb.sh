#!/bin/bash

ligands=( 1 2 3 4 5 6 7 9 10 11 12 13 14 15 )

parent_dir="$HOME/projects/metalloenzymes/nonbonded_model_vim2/"
prep_dir="$parent_dir/system_preparation"
steps=( 1 2 )

for i in ${ligands[@]}
do
	ligand_directory="$parent_dir/ligand_$i"
	cd $ligand_directory

	step_4n1_directory="$ligand_directory/step_4n1"
	cd $step_4n1_directory
	echo $step_4n1_directory
	sed -i "s/L8J/LIG/g" vim2.in
	sed -i "s/vim2.pdb/prepared_system.pdb/g" vim2.in
	# MCPB.py -i vim2.in -s 1m 
	# MCPB.py -i vim2.in -s 4n1

	step_4n2_directory="$ligand_directory/step_4n2"
	cd $step_4n2_directory
	echo $step_4n2_directory
	sed -i "s/L8J/LIG/g" vim2.in
	sed -i "s/vim2.pdb/prepared_system.pdb/g" vim2.in
	MCPB.py -i vim2.in -s 4n2
	sed -i "2 isource leaprc.gaff2" vim2_tleap.in
	sed -i "s/LIG.mol2/LIG/g" vim2_tleap.in
	sed -i "s/WAT.mol2/WAT/g" vim2_tleap.in

done

