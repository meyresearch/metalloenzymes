#!/bin/bash
#SBATCH --job-name production
#SBATCH --output equil.%j.out
#SBATCH --error equil.%j.err
#SBATCH --mem 4096


$AMBERHOME/bin/pmemd.cuda -O -i md.in -p vim2_solv.prmtop -c 9md.rst7\
 -ref 9md.rst7 -o vim2_solv_md.mdout -r vim2_solv_md.rst7 -x vim2_solv_md.nc
