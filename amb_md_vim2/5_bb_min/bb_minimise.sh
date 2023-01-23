#!/bin/bash

export CUDA_VISIBLE_DEVICES=0

pmemd.cuda -O -i bb_min.in -o bb_min.out -p ../1_min/vim2_solv.prmtop -c ../4_lower/lower.rst7 -r bb_min.rst7 -inf bb_min.info -ref ../4_lower/lower.rst7 -x mdcrd.bb_min
