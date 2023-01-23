#!/bin/bash

ligands=( 1 2 3 4 5 6 7 9 10 11 12 13 14 15 16 )

cd ..
for i in ${ligands[@]}
do
	cd ligand_$i
	echo "ligand_$i"
	rsync -av azuma:/home/jguven/projects/metalloenzymes/dft/ligand_$i/* .
	rsync -av eddie:/exports/eddie/scratch/s1607171/ligand_$i/* .
	cd ..
done

