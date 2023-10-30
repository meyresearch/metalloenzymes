#!/bin/bash
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=1
#SBATCH --mem=4GB

heat_dir=HEAT_DIRECTORY

mpirun -v sander.MPI -O -i min.cfg -o min.out -p min.prm7 -c min.rst7 -r min.rst7 -inf min.info
cp min.rst7 $heat_dir/heat.rst7
