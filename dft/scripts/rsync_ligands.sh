#!/bin/bash

cd ..
ligands=( 1 2 3 4 5 6 7 9 10 11 12 13 14 15 16 )

for i in ${ligands[@]}
do
	rsync -av ligand_$i/ligand_$i.pdb azuma:projects/metalloenzymes/dft/ligand_$i/ 
done
