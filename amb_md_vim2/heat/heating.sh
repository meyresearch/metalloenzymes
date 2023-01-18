#!/bin/bash

export CUDA_VISIBLE_DEVICES=0

pmemd.cuda -O -i heat.in -o heat.out -p ../min/vim2_solv.prmtop -c ../min/min.rst7 -r heat.rst7 -inf heat.info -ref ../min/min.rst7 -x mdcrd.heat 
