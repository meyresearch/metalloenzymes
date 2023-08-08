#!/bin/bash

ligands=( 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 )

parent_dir="$HOME/projects/alchemistry/vim2_partially_protonated_ligand/parameterisation"
scripts_dir="$HOME/projects/metalloenzymes/scripts/"
for i in {1..16}
do
	ligand_directory="$parent_dir/ligand_$i"
	cd $ligand_directory
	cp $scripts_dir/small_opt.sh .
	sed -i "s/#$ -N gssn16/#$ -N sopt_$i/g" small_opt.sh
	sed -i "s/%Mem=3000MB/%Mem=12000MB/g" vim2_small_opt.com
	sed -i "s/%NProcShared=2/%NProcShared=8/g" vim2_small_opt.com
	rsync -av small_opt.sh eddie:/exports/eddie/scratch/s1607171/part_prot/ligand_$i/
	rsync -av vim2_small_opt.com eddie:/exports/eddie/scratch/s1607171/part_prot/ligand_$i/
	
done
