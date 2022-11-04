#!/bin/bash 
cd ../runs/equilibration/unbound/ligand_10/r_npt 
gmx grompp -f r_npt.mdp -c ../nvt/nvt.gro -r ../nvt/nvt.gro -p r_npt.top -t ../nvt/nvt.cpt -o r_npt.tpr 
nohup gmx mdrun -v -deffnm r_npt -nt 1 -nb gpu > nohup_r_npt.out &