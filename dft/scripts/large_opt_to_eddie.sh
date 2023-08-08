#!/bin/bash

ligands=( 2 3 4 6 7 9 10 14 15 16 )

cd ..
for i in ${ligands[@]}
do
	cd ligand_$i/
	cp ../scripts/large_opt.sh .
	sed -i "s/#$ -N large_opt/#$ -N lopt_$i/g" large_opt.sh
	cp vim2_large_mk.com vim2_large_opt.com
	sed -i "s/%Mem=3000MB/%Mem=12000MB/g" vim2_large_opt.com
	sed -i "s/%NProcShared=2/%NProcShared=8/g" vim2_large_opt.com
	sed -i "s/ Pop(MK,ReadRadii)/ /g" vim2_large_opt.com
	rsync -av vim2_large_opt.com eddie:/exports/eddie/scratch/s1607171/ligand_$i/
	rsync -av large_opt.sh eddie:/exports/eddie/scratch/s1607171/ligand_$i/
	cd ../
done
