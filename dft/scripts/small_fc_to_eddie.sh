#!/bin/bash

ligands_1=( 1 4 5 6 7 9 10 11 13 14 15 )

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

ligands_2=( 2 3 16 )

for i in ${ligands_2[@]}
do
	echo "ligand $i"
	cd ligand_$i/
	cp ../scripts/small_fc.sh .
	sed -i "s/#$ -N gssn16/#$ -N sfc_$i/g" small_fc.sh
	sed -i "s/%Mem=3000MB/%Mem=12000MB/g" vim2_small_fc.com
	sed -i "s/%NProcShared=2/%NProcShared=8/g" vim2_small_fc.com
	sed -i "s/%Chk=vim2_small_opt_2.chk/%Chk=vim2_small_opt.chk/g" vim2_small_fc.com
	rsync -av vim2_small_fc.com eddie:/exports/eddie/scratch/s1607171/ligand_$i/
	rsync -av small_fc.sh eddie:/exports/eddie/scratch/s1607171/ligand_$i/
	cd ../
done
