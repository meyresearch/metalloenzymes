#!/bin/bash

ligands=( 1 2 3 4 5 6 7 9 10 11 12 13 14 15 16 )

cd ..
for i in ${ligands[@]}
do
	cd ligand_$i/
	echo "====================ligand $i==================="
	# Step 2 - Seminario method
        MCPB.py -i vim2.in -s 2
	# Step 3 - RESP charge fitting
	MCPB.py -i vim2.in -s 3
	# Step 4 - Generate tleap input
	MCPB.py -i vim2.in -s 4
	# Edit vim2_tleap.in	
	sed -i '1s/^/source leaprc.gaff2\n/' vim2_tleap.in
	# Run tleap
	tleap -s -f vim2_tleap.in > vim2_tleap.out
	cd ..
	echo "================================================"
done
