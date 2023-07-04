#!/bin/bash

network_file=$1

# Open network.dat file and get data out
transformations=()
lambdas=()
engines=()

while IFS= read -r lines
do
	while read -a line
		do
			transformation=${line[0]}~${line[1]}
			n_windows=${line[2]}
			engine=${line[-1]}
			transformations+=('$transformation')
			lambdas+=('$n_windows')
			engines+=('$engine')
			echo "trans: $transformation, l: $n_windows, engine: $engine"
		done <<< $lines
done < $network_file
