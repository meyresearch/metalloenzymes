#!/bin/bash

ligands=( 1 2 3 4 5 6 7 9 10 11 12 13 14 15 )

parent_dir="$HOME/projects/metalloenzymes/nonbonded_model_vim2/"

for i in ${ligands[@]}
do
	ligand_directory="$parent_dir/ligand_$i"
	production_dir="$ligand_directory/step_4n2/production"
	rm -r $production_dir
done