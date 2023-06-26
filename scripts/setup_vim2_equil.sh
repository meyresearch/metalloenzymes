#!/bin/bash

ligands=( 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 )

afe_dir="/home/jguven/projects/metalloenzymes/kpc2/afe/"
# network="$afe_dir/network_fwd.dat"
afe_dir="$HOME/projects/metalloenzymes/vim_2/inputs/"
parent_dir="$HOME/projects/metalloenzymes/nonbonded_model_vim2/"
if [[ ! -d $afe_dir ]]; then
    mkdir $afe_dir
    mkdir "$afe_dir/ligands/"
    mkdir "$afe_dir/protein/"
fi

for i in ${ligands[@]}
do
    old_ligand_directory="$parent_dir/ligand_$i"
    ligand_file="$old_ligand_directory/L8J.pdb"
    prmtop="$old_ligand_directory/vim2_solv.prmtop"
    rst7="$old_ligand_directory/production/vim2_md.rst7"

    cp $ligand_file $afe_dir/ligands/ligand_${i}.pdb
    cp $prmtop $afe_dir/protein/system_${i}.prm7
    cp $rst7 $afe_dir/protein/system_${i}.rst7

done