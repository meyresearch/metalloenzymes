#!/bin/bash
ligands=($(seq 1 1 16))
for i in ${ligands[@]}; 
do 
	cp ../../first_year_project/kpc2/ligands/ligand_$i/lig_h_$i.mol2 ../all_ligands/ligand_$i.mol2
done
