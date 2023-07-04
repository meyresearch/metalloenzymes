#!/bin/bash

ligands=( 1 3 4 6 7 9 10 12 13 16 )
dft_dir=$HOME/projects/metalloenzymes/dft/
for i in ${ligands[@]}
do
	cd $dft_dir/ligand_$i
	if [[ ! -d $dft_dir/ligand_$i/gromacs ]]; then
		mkdir $dft_dir/ligand_$i/gromacs
	fi
	echo $PWD
	mv 01_minimisation $dft_dir/ligand_$i/gromacs/01_minimisation
	mv 02_nvt $dft_dir/ligand_$i/gromacs/02_nvt
	mv 03_npt $dft_dir/ligand_$i/gromacs/03_npt
	mv 04_md $dft_dir/ligand_$i/gromacs/04_md
done
