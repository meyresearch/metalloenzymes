#!/bin/bash 
cd ../runs/equilibration/unbound/ligand_8/r_nvt 
gmx grompp -f r_nvt.mdp -c ../min/min.gro -r ../min/min.gro -p r_nvt.top -o r_nvt.tpr 
nohup gmx mdrun -v -deffnm r_nvt -nt 1 -nb gpu > nohup_r_nvt.out &