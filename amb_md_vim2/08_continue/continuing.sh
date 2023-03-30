#!/bin/bash

export CUDA_VISIBLE_DEVICES=0

pmemd.cuda -O -i continue.in -o continue.out -p ../1_min/vim2_solv.prmtop -c ../7_reduce/reduce.rst7 -r continue.rst7 -inf continue.info -ref ../7_reduce/reduce.rst7 -x mdcrd.continue
