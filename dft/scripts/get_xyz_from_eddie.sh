#!/bin/bash

ligands=( 1 2 3 4 5 6 7 9 10 11 12 13 14 15 16 )
ligands=( 2 5 11 14 15 )

for i in ${ligands[@]}
do
	rsync -av  eddie:/exports/eddie/scratch/s1607171/ligand_$i/small_ligand_$i.xyz ../geometries/
done
