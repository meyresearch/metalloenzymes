#!/bin/bash

ligands=( 2 5 11 14 15 )

cd ..
for i in ${ligands[@]}
do
	cd ligand_$i/
	cp ../scripts/small_opt.sh .
	sed -i "s/#$ -N gssn16/#$ -N sopt_$i/g" small_opt.sh
	sed -i "s/%Mem=3000MB/%Mem=12000MB/g" vim2_small_opt.com
	sed -i "s/%NProcShared=2/%NProcShared=8/g" vim2_small_opt.com
	rsync -av * eddie:/exports/eddie/scratch/s1607171/ligand_$i/
	cd ../
done
