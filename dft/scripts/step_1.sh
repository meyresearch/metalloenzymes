#!/bin/bash

#ligands=( 1 2 3 4 5 6 7 9 10 11 12 13 14 15 16 )
ligands=( 2 5 11 14 15 )
cd .. 
for i in ${ligands[@]}
do
	cp protein/vim2.in ligand_$i/
	cd ligand_$i/
	MCPB.py -i vim2.in -s 1
	cd ..
done
