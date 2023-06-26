#!/bin/bash
#SBATCH -o /home/jguven/projects/metalloenzymes/slurm_logs/plot_output.%A.%a.slurm.out
#SBATCH -e /home/jguven/projects/metalloenzymes/slurm_logs/plot_output.%A.%a.slurm.err
# Only launch a max of 1 task
#SBATCH -n 1
# allocate 1 gpu per job
#SBATCH --gres=gpu:1
# allocate 8 processors per gpu
#SBATCH --cpus-per-gpu 10
#SBATCH --mem 4096
# set job name
#SBATCH --job-name=plot
date

idx=$SLURM_ARRAY_TASK_ID

python analyse_all_md.py --ligand-number $idx
