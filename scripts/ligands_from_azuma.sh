#!/bin/bash

ligands=($(seq 2 1 16))
for i in ${ligands[@]};
do
        rsync azuma:projects/first_year_project/kpc2/ligands/ligand_$i/lig_h_$i.mol2 ../all_ligands/ligand_$i.mol2
done
