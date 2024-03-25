#! /bin/bash

dos2unix "$transformations_file"


transformations_array=()
lambdas_array=()
n_windows_array=()

while IFS="," read -a line
do
    transformation=${line[1]}~${line[2]}
    lambda=${line[5]}
    n_windows=${line[4]}
    transformations_array+=("$transformation")
    lambdas_array+=("$lambda")
    n_windows_array+=("$n_windows")
done < <(tail -n +2 "$transformations_file")
