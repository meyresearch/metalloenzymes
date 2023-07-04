#!/bin/bash

ligands=( 1 2 3 4 5 6 7 9 10 11 12 13 14 15 16 )

for i in ${ligands[@]}
do
	rm ligand_$i/vim2_with_$i.pdb
	#mv ligand_$i/*.pdb ligand_$i/vim2_with_$i.pdb
	#echo "mv ligand_$i/*.pbd ligand_$i/vim_2_with_$i.pdb"
done
