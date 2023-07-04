#!/bin/bash

for i in range $(seq 1 1 16)
do
	rsync ../kpc2/inputs/ligands/docked_$i.sdf pluto:projects/metalloenzymes/kpc2/inputs/ligands/
done
