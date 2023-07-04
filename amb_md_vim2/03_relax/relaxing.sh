#!/bin/bash

export CUDA_VISIBLE_DEVICES=0

pmemd.cuda -O -i relax.in -o relax.out -p ../1_min/vim2_solv.prmtop -c ../2_heat/heat.rst7 -r relax.rst7 -inf relax.info -ref ../2_heat/heat.rst7 -x mdcrd.relax
