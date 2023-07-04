#!/bin/bash

ligands=( 1 3 4 6 7 9 10 12 13 16 )

cd ..
for i in ${ligands[@]}
do
	echo "ligand $i"
	cd ligand_$i/
	amb2gro_top_gro.py -p vim2_solv.prmtop -c vim2_solv.inpcrd -t vim2_solv.top -g vim2_solv.gro -b vim2_solv.gromacs.pdb
	cd ..
done
