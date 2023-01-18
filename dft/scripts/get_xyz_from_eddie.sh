#!/bin/bash

ligands=( 1 2 3 4 5 6 7 9 10 11 12 13 14 15 16 )


for i in ${ligands[@]}
do
	rsync -av  eddie:/exports/eddie/scratch/s1607171/ligand_$i/ligand_$i.xyz ../geometries/
done
