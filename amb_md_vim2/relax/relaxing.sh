#!/bin/bash

export CUDA_VISIBLE_DEVICES=0

pmemd.cuda -O -i relax.in -o relax.out -p ../min/vim2_solv.prmtop -c ../heat/heat.rst7 -r relax.rst7 -inf relax.info -ref ../heat/heat.rst7 -x mdcrd.relax
