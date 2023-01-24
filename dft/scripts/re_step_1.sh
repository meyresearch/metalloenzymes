#!/bin/bash

ligands=( 2 5 11 14 15 )

cd ..
for i in ${ligands[@]}
do
	cp protein/vim2.in ligand_$i/
	cp protein/ZN1.mol2 ligand_$i/
	cp protein/ZN2.mol2 ligand_$i/
	cp protein/WAT.mol2 ligand_$i/
	cd ligand_$i/
	MCPB.py -i vim2.in -s 1
	cd ..
done
