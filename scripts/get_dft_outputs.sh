#!/bin/bash

ligands=( 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 )
parent_dir="$HOME/projects/alchemistry/vim2_deprotonated_ligand/parameterisation/"

for i in ${ligands[@]}
do
	lig_dir=$parent_dir/ligand_$i
	cd $lig_dir
	rsync -av eddie:/exports/eddie/scratch/s1607171/ligand_$i/vim2_small_opt.log .
	rsync -av eddie:/exports/eddie/scratch/s1607171/ligand_$i/vim2_small_opt.chk .
done

