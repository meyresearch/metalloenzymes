#!/bin/bash

ligands=( 1 2 3 4 5 6 7 9 10 11 12 13 14 15 )

parent_dir="$HOME/projects/metalloenzymes/nonbonded_model_vim2/"

for i in ${ligands[@]}
do
	ligand_directory="$parent_dir/ligand_$i"
	step_4n2_directory="$ligand_directory/step_4n2"
	cd $step_4n2_directory
	echo $step_4n2_directory
    tleap -s -f vim2_tleap.in > vim2_tleap.out
done