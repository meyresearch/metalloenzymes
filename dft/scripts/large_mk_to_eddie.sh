#!/bin/bash

ligands=( 1 2 3 4 5 6 7 9 10 11 12 13 14 15 16 )
ligands=( 2 5 11 14 15 )
cd ..
for i in ${ligands[@]}
do
	cd ligand_$i/
	cp ../scripts/large_mk.sh .
	sed -i "s/#$ -N gssn16/#$ -N mk_$i/g" large_mk.sh
	sed -i "s/%Mem=3000MB/%Mem=12000MB/g" vim2_large_mk.com
	sed -i "s/%NProcShared=2/%NProcShared=8/g" vim2_large_mk.com
	rsync -av vim2_large_mk.com eddie:/exports/eddie/scratch/s1607171/ligand_$i/
	rsync -av large_mk.sh eddie:/exports/eddie/scratch/s1607171/ligand_$i/
	cd ../
done
