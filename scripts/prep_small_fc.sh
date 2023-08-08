#!/bin/bash


parent_dir="$HOME/projects/alchemistry/vim2_partially_protonated_ligand/parameterisation"
scripts_dir="$HOME/projects/metalloenzymes/scripts/"
for i in {9..16}
do
	ligand_directory="$parent_dir/ligand_$i"
	echo $ligand_directory
	cd $ligand_directory
	cp $scripts_dir/small_fc.sh .
	sed -i "s/#$ -N sfc_3/#$ -N sfc_$i/g" small_fc.sh
	sed -i "s/%Mem=3000MB/%Mem=12000MB/g" vim2_small_fc.com
	sed -i "s/%NProcShared=2/%NProcShared=8/g" vim2_small_fc.com
	rsync -av small_fc.sh eddie:/exports/eddie/scratch/s1607171/part_prot/ligand_$i/
	rsync -av vim2_small_fc.com eddie:/exports/eddie/scratch/s1607171/part_prot/ligand_$i/
	
done
