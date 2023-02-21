#!/bin/bash
ligands=($(seq 1 1 16))
for i in ${ligands[@]};
do
        rsync -av azuma:projects/first_year_project/kpc2/ligands/ligand_$i/lig_h_$i.mol2 $HOME/projects/metalloenzymes/all_ligands/vim2/ligand_$i.mol2
done

