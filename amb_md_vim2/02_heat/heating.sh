#!/bin/bash

export CUDA_VISIBLE_DEVICES=0

pmemd.cuda -O -i heat.in -o heat.out -p ../1_min/vim2_solv.prmtop -c ../1_min/min.rst7 -r heat.rst7 -inf heat.info -ref ../1_min/min.rst7 -x mdcrd.heat 
