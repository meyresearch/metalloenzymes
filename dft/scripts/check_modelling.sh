#!/bin/bash

ligands=( 1 2 3 4 5 6 7 9 10 11 12 13 14 15 16 )

cd ..
for i in ${ligands[@]}
do
	cd ligand_$i/
	cp ../scripts/parmed_*.in .
	echo "====================ligand $i===================="
	cpptraj -p vim2_solv.prmtop > cpptraj.out &
	parmed -i parmed_1.in -p vim2_solv.prmtop > parmed_1.out &
	parmed -i parmed_2.in -p vim2_solv.prmtop > parmed_2.out &
	cd ..
	echo "================================================"
done
