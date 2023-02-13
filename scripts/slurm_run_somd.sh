#!/bin/bash
# sys arg: $1: network file
#afe_dir=$HOME/projects/metalloenzymes/kpc2/afe


network=$1
mapfile network_file < $network
for perturbation in "${network_file[@]}"
do
	IFS=' ' read -a ligpair <<< $perturbation
	n_windows=$(( ${ligpair[2]} ))
	echo "--Processing ${ligpair[0]} ${ligpair[1]} with $n_windows windows---" 	
	# ligpair[0], ligpair[1]: ligand 1, ligand 2
	# ligpair[2]: number of lambda windows
	# ligpair[3]: lambda array
	# ligpair[4]: engine
	sbatch --array=0-$(($n_windows-1)) somd_afe.sh ${ligpair[0]} ${ligpair[1]} ${ligpair[4]} ${ligpair[3]}
done