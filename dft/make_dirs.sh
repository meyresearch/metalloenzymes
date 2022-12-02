#!/bin/bash

ligands=( 1 2 3 4 5 6 7 9 10 11 12 13 14 15 16 )
for i in ${ligands[@]}
do
	mkdir -p ligand_$i
	touch ligand_$i/.placeholder
done

