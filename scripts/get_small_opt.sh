#!/bin/bash

for i in {9..16}
do
	rsync -av eddie:/exports/eddie/scratch/s1607171/part_prot/ligand_${i}/lig_${i}.xyz ~/projects/alchemistry/vim2_partially_protonated_ligand/parameterisation/ligand_${i}/
	rsync -av eddie:/exports/eddie/scratch/s1607171/part_prot/ligand_$i/vim2_small_opt.* ~/projects/alchemistry/vim2_partially_protonated_ligand/parameterisation/ligand_${i}/
	#mv ~/projects/qmmm/nonbonded_models/integer_charge/part_protonated/ligand_${i}/vim2_small_opt.chk /home/jguven/projects/alchemistry/vim2_partially_protonated_ligand/parameterisation/ligand_${i}/
	#mv ~/projects/qmmm/nonbonded_models/integer_charge/part_protonated/ligand_${i}/vim2_small_opt.log /home/jguven/projects/alchemistry/vim2_partially_protonated_ligand/parameterisation/ligand_${i}/
	#mv ~/projects/qmmm/nonbonded_models/integer_charge/part_protonated/ligand_${i}/vim2_small_opt.com /home/jguven/projects/alchemistry/vim2_partially_protonated_ligand/parameterisation/ligand_${i}/
	#mv ~/projects/qmmm/nonbonded_models/integer_charge/part_protonated/ligand_${i}/lig_${i}.xyz /home/jguven/projects/alchemistry/vim2_partially_protonated_ligand/parameterisation/ligand_${i}/

done
