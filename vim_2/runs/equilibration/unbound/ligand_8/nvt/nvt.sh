#!/bin/bash 
cd ../runs/equilibration/unbound/ligand_8/nvt 
gmx grompp -f nvt.mdp -c ../r_nvt/r_nvt.gro -p nvt.top -t ../r_nvt/r_nvt.cpt -o nvt.tpr 
nohup gmx mdrun -v -deffnm nvt -nt 1 -nb gpu > nohup_nvt.out &