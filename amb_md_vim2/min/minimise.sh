#!/bin/bash

export CUDA_VISIBLE_DEVICES=0

# -O overwrite output files, -i, -o input and output files
# -p, -c prm7, rst7
# -r restart output file
# -inf write MD info file with simulation status
# -ref reference coordinates: same as input coordinates
# -x file with trajectory

pmemd.cuda -O -i min.in -o min.out -p vim2_solv.prmtop -c vim2_solv.inpcrd -r min.rst7 -inf min.info -ref vim2_solv.inpcrd -x mdcrd.min 
