#!/bin/bash

ligands=(  1 2 3 4 5 6 7 9 10 11 12 13 14 15 )

parent_dir="$HOME/projects/metalloenzymes/nonbonded_model_vim2/"
# prep_dir="$parent_dir/system_preparation"
# steps=( 1 2 )

for i in ${ligands[@]}
do
	ligand_directory="$parent_dir/ligand_$i"
	cd $ligand_directory
	# for j in ${steps[@]}
	# do
		# step_directory="$ligand_directory/step_4n$j"
		# cd $step_directory
		# echo $step_directory

	cat vim2_protonated.pdb  L8J.pdb HOH.pdb > vim2_cat.pdb
	pdb4amber vim2_cat.pdb > vim2_fixed.pdb
		
		
	# done
done

