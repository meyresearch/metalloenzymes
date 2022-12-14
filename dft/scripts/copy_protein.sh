#!/bin/bash

cd ..

ligands=( 1 2 3 4 5 6 7 9 10 11 12 13 14 15 16 )

for i in ${ligands[@]}
do 
	cp protein/vim2_split.pdb ligand_$i/
        cp protein/WAT.tleap.pdb ligand_$i/
	cp protein/*.mol2 ligand_$i/
done	
