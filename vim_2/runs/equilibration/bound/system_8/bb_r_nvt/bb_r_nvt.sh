#!/bin/bash 
cd ../runs/equilibration/bound/system_8/bb_r_nvt 
gmx grompp -f bb_r_nvt.mdp -c ../min/min.gro -r ../min/min.gro -p bb_r_nvt.top -o bb_r_nvt.tpr 
nohup gmx mdrun -v -deffnm bb_r_nvt -nt 1 -nb gpu > nohup_bb_r_nvt.out &