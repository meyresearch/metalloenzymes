#!/bin/bash
#$ -N nvt

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
gmx grompp -f nvt.mdp -c em.gro -r em.gro -p test.top -o nvt.tpr
gmx mdrun -deffnm nvt -nb gpu
