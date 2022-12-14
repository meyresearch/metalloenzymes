#!/bin/bash

cd ..

ligands=( 1 2 3 4 5 6 7 9 10 11 12 13 14 15 16 )

for i in ${ligands[@]}
do 
	cd ligand_$i
	cat vim2_split.pdb LIG.pdb WAT.tleap.pdb > vim2_complex.pdb
	pdb4amber -i vim2_complex.pdb -o vim2_complex.amber.pdb
	cd ..
done	
