#!/bin/bash

export CUDA_VISIBLE_DEVICES=0

pmemd.cuda -O -i bb_relax.in -o bb_relax.out -p ../1_min/vim2_solv.prmtop -c ../5_bb_min/bb_min.rst7 -r bb_relax.rst7 -inf bb_relax.info -ref ../5_bb_min/bb_min.rst7 -x mdcrd.bb_relax
