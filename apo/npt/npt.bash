#!/bin/bash
#$ -N npt

# Working directory
#$ -cwd

# Runtime limit
#$ -l h_rt=48:00:00

# Request system ram
#$ -l h_vmem=12G

# Request GPUs
#$ -pe gpu-titanx 2

# Initialise environmental modules. Load gromacs
. /etc/profile.d/modules.sh
module load chem/gromacs/2021.2

# Commands
gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -p test.top -o npt.tpr
gmx mdrun -deffnm npt -nb gpu
