#! /bin/bash

dos2unix "$ligands_dat_file"
dos2unix "$transformations_file"

ligand_array=()
while read line
do 
    clean=$(sed "s/\r$//" <<< $ligands)
    ligand_array+=("$line")

done < $ligands_dat_file

transformations_array=()
lambdas_array=()
n_windows_array=()

while IFS="," read -a line
do
    transformation=${line[1]}~${line[3]}
    lambda=${line[7]}
    n_windows=${line[6]}
    transformations_array+=("$transformation")
    lambdas_array+=("$lambda")
    n_windows_array+=("$n_windows")
done < <(tail -n +2 "$transformations_file")
