#!/bin/bash
#SBATCH --job-name=geo_opt
#SBATCH --mem 4096
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30

srun python geo_opt.py
