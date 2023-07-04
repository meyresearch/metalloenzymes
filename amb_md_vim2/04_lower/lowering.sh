#!/bin/bash

export CUDA_VISIBLE_DEVICES=0

pmemd.cuda -O -i lower.in -o lower.out -p ../1_min/vim2_solv.prmtop -c ../3_relax/relax.rst7 -r lower.rst7 -inf lower.info -ref ../3_relax/relax.rst7 -x mdcrd.lower
