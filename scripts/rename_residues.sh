#!/bin/bash

parent_dir=$1

for i in {1..16}
do
	ligand_directory="$parent_dir/ligand_$i"
	cd $ligand_directory

	sed -i "s/HETATM/ATOM  /g" vim2_fixed.pdb
	sed -i "s/TYN/TYR/g" vim2_fixed.pdb
	sed -i "s/ARN/ARG/g" vim2_fixed.pdb
		
done

