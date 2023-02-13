#!/bin/bash

ligands_1=( 2 5 11 14 15 )

cd ..
for i in ${ligands_1[@]}
do
	echo "ligand $i"
	cd ligand_$i/
	cp ../scripts/small_fc.sh .
	sed -i "s/#$ -N gssn16/#$ -N sfc_$i/g" small_fc.sh
	sed -i "s/%Mem=3000MB/%Mem=12000MB/g" vim2_small_fc.com
	sed -i "s/%NProcShared=2/%NProcShared=8/g" vim2_small_fc.com
	rsync -av vim2_small_fc.com eddie:/exports/eddie/scratch/s1607171/ligand_$i/
	rsync -av small_fc.sh eddie:/exports/eddie/scratch/s1607171/ligand_$i/
	cd ../
done
