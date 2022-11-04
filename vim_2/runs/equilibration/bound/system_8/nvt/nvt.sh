#!/bin/bash 
cd ../runs/equilibration/bound/system_8/nvt 
gmx grompp -f nvt.mdp -c ../bb_r_nvt/bb_r_nvt.gro -p nvt.top -t ../bb_r_nvt/bb_r_nvt.cpt -o nvt.tpr 
nohup gmx mdrun -v -deffnm nvt -nt 1 -nb gpu > nohup_nvt.out &