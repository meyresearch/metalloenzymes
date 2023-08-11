#!/bin/bash
#SBATCH --output logs/meze.%a.out
#SBATCH --error logs/meze.%a.err
#SBATCH -n 1
#SBATCH --cpus-per-task 10 
#SBATCH --mem 4096
#SBATCH --job-name=meze


idx=$SLURM_ARRAY_TASK_ID

python ~/projects/metalloenzymes/meze/meze.py --step 2 --group-name kpc2 --input-pdb-file inputs/protein/kpc2.input.pdb --ligand-charge -1 --ligand-index $idx  

