#!/bin/bash

gpu=$1

export CUDA_VISIBLE_DEVICES=$gpu

pmemd.cuda -O -i md.in -o md.out -p vim2_equil.prmtop -c vim2_equil.rst7 -ref vim2_equil.rst7 -r vim2_equil.rst7 -x md.nc 
