#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=28
#SBATCH --mem=10GB
#SBATCH --time=05:00:00

mpirun -np 28 sander.MPI -O -i min-qmmm.in -o min-qmmm.out -c vim2_solv.inpcrd -r vim2_solv.inpcrd -p vim2_solv.prmtop

