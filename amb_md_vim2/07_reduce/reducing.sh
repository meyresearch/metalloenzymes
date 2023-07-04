#!/bin/bash

export CUDA_VISIBLE_DEVICES=0

pmemd.cuda -O -i reduce.in -o reduce.out -p ../1_min/vim2_solv.prmtop -c ../6_bb_relax/bb_relax.rst7 -r reduce.rst7 -inf reduce.info -ref ../6_bb_relax/bb_relax.rst7 -x mdcrd.reduce
