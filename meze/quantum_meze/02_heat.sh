#!/bin/bash
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=1
#SBATCH --mem=4GB

outputs_directory=OUTPUTS_DIRECTORY

mpirun -v sander.MPI -O -i heat.cfg -o heat.out -p heat.prm7 -c heat.rst7 -r heat.rst7 -inf heat.info
cp heat.rst7 $outputs_directory/qmmm.rst7
