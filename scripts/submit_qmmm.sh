#!/bin/bash

export AMBERHOME=/home/jguven/Software/amber22/

parent_dir=$HOME/projects/qmmm/nonbonded_models/integer_charge/part_protonated/
qmmm_dir=$HOME/projects/qmmm/qm_input_files/

for i in {1..16}
do
	echo ligand_$i
	ligand_dir=$parent_dir/ligand_$i
	cd $ligand_dir
	cp $qmmm_dir/minimise.sh .
	cp $qmmm_dir/equilibrate.sh .
	cp $qmmm_dir/production.sh .
	cp $qmmm_dir/qmmm.sh .
	cp $qmmm_dir/min.in .
	cp $qmmm_dir/equil.in .
	cp $qmmm_dir/prod.in .
	sbatch qmmm.sh
	sleep 10
done

