#!/bin/bash

# -o controls output. %A is replaced with job ID and %a with array index
#SBATCH -o /home/jguven/projects/metalloenzymes/slurm_logs/equil.%a.slurm.out
#SBATCH -e /home/jguven/projects/metalloenzymes/slurm_logs/equil.%a.slurm.err
# Only launch a max of 1 task
#SBATCH -n 1
# allocate 1 gpu per job
#SBATCH --gres=gpu:1
# allocate 8 processors per gpu
#SBATCH --cpus-per-gpu 10
#SBATCH --mem 4096
# set job name
#SBATCH --job-name=equilibration
# specify number of threads
export OMP_NUM_THREADS=8

ligand_file=$1

readarray -t ligands < $ligand_file

ligand=${ligands[$SLURM_ARRAY_TASK_ID]}

echo $ligand

slurm_logs_dir=/home/jguven/projects/metalloenzymes/slurm_logs/
if [[ ! -d $slurm_logs_dir ]]; then
       mkdir $slurm_logs_dir
fi

python mequilibrafe.py /home/jguven/projects/metalloenzymes/vim_2/ $SLURM_ARRAY_TASK_ID
#srun python mequilibrafe.py kpc2 ../azetidinimines/kpc2/ $ligand 
