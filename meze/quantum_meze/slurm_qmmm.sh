#!/bin/bash

min_dir=MINIMISATION_DIRECTORY
heat_dir=HEAT_DIRECTORY
outputs_dir=OUTPUTS_DIRECTORY



cd min
min=$(sbatch minimise.sh)
cd $HOME/projects/alchemistry/vim2/deprotonated_ligand/qmmm/equilibration/ligand_1/heat/
heat=$(sbatch --dependency=afterok:${min} equilibrate.sh)

