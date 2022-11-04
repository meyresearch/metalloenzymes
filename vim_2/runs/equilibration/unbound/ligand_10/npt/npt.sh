#!/bin/bash 
cd ../runs/equilibration/unbound/ligand_10/npt 
gmx grompp -f npt.mdp -c ../r_npt/r_npt.gro -p npt.top -t ../r_npt/r_npt.cpt -o npt.tpr 
nohup gmx mdrun -v -deffnm npt -nt 1 -nb gpu > nohup_npt.out &