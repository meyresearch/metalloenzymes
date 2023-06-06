#!/bin/bash

ligands=( 1 2 3 4 5 6 7 9 10 11 12 13 14 15 )

parent_dir="$HOME/projects/metalloenzymes/nonbonded_model_vim2/"
prep_dir="$parent_dir/system_preparation"
steps=( 1 2 )

for i in ${ligands[@]}
do
	ligand_directory="$parent_dir/ligand_$i"
	cd $ligand_directory
	for j in ${steps[@]}
	do
		step_directory="$ligand_directory/step_4n$j"
		cd $step_directory
		echo $step_directory
		cp $prep_dir/ZN* .
		grep -v '^HETATM' vim2.pdb > protein.pdb
		cat protein.pdb ZN1.pdb ZN2.pdb LIG.pdb WAT.pdb > system.pdb
		pdb4amber system.pdb > prepared_system.pdb
		sed -i "s/HOH/WAT/g" prepared_system.pdb
		sed -i "s/HOH/WAT/g" WAT.mol2
	done
done

