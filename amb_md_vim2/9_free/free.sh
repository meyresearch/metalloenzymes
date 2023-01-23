#!/bin/bash

export CUDA_VISIBLE_DEVICES=0

pmemd.cuda -O -i free.in -o free.out -p ../1_min/vim2_solv.prmtop -c ../8_continue/continue.rst7 -r free.rst7 -inf free.info -ref ../8_continue/continue.rst7 -x mdcrd.free
