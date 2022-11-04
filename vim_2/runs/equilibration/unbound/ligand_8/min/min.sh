#!/bin/bash 
cd ../runs/equilibration/unbound/ligand_8/min 
gmx grompp -f min.mdp -c min.gro -p min.top -o min.tpr 
nohup gmx mdrun -v -deffnm min -nt 1 -nb gpu > nohup_min.out &