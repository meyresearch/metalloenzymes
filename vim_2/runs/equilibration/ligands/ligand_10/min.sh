#!/bin/bash 
cd ../runs/equilibration/ligands/ligand_10 
gmx grompp -f min.mdp -c min.gro -p min.top -o min.tpr 
nohup gmx mdrun -v -deffnm min -nt 1 -nb gpu > nohup_min.out &