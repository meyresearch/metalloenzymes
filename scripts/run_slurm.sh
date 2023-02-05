#!/bin/bash

# Make directory for slurm out and err files
if [[ ! -d ../slurm_logs ]]; then
	mkdir ../slurm_logs
fi

network_file=$1

# Open network.dat file and get data out
transformations=()
lambdas=()
engines=()

sed -i "s/,//g" $network_file
while IFS= read -r lines
do
	while read -a line
		do
			transformation=${line[0]}~${line[1]}
			n_windows=${line[2]}
			engine=${line[-1]}
			transformations+=("$transformation")
			lambdas+=("$n_windows")
			engines+=("$engine")
#			echo "trans: $transformation, l: $n_windows, engine: $engine"
		done <<< $lines
done < $network_file

echo ${transformations[@]}
echo ${lambdas[@]}
echo ${engines[@]}

for i in "${!transformations[@]}"
do
	#jid_afe=$(sbatch --parsable --array=0-$((${lambdas[i]}-1)) slurm_run_somd.sh ${transformations[i]} ${lambdas[i]})
	echo "afe for ${transformations[i]} with ${lambdas[i]} windows"
done
