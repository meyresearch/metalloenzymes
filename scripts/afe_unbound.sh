#!/bin/bash

gpu=$1
unbound_data=$2

echo "Running on GPU $gpu"

while IFS= read -r unbound_transformation
do
	cd $unbound_transformation
	# MINIMISATION
	echo "Starting minimisation..."
       	cd min

	gmx grompp -f min.mdp -c min.gro -p min.top -o min.tpr
	gmx mdrun -v -deffnm min -gpu_id $gpu -c ../r_nvt/r_nvt.gro
	
	sleep 10
	cp min.top ../r_nvt/r_nvt.top
	cd ../
	echo "Minimisation complete."
	
	# Restrained NVT
	echo "Starting restrained NVT equilibration..."
	cd r_nvt
	gmx grompp -f r_nvt.mdp -c r_nvt.gro -r r_nvt.gro -p r_nvt.top -o r_nvt.tpr -maxwarn 1
	gmx mdrun -v -deffnm r_nvt -gpu_id $gpu -c ../nvt/nvt.gro
	
	sleep 10
	cp r_nvt.top ../nvt/nvt.top
	cd ../
	echo "Restrained NVT complete."

	# NVT
	echo "Starting NVT equilibration..."
	cd nvt
	gmx grompp -f nvt.mdp -c nvt.gro -t ../r_nvt/r_nvt.cpt -p nvt.top -o nvt.tpr -maxwarn 1
	gmx mdrun -v -deffnm nvt -gpu_id $gpu -c ../npt/npt.gro

	sleep 10
	cp nvt.top ../npt/npt.top
	cd ../
	echo "NVT equilibration complete."

	# NPT
	echo "Starting NPT equilibration..."
	cd npt
	gmx grompp -f npt.mdp -c npt.gro -p npt.top -t ../nvt/nvt.cpt -o npt.tpr -maxwarn 1
	gmx mdrun -v -deffnm npt -gpu_id $gpu -c ../afe/gromacs.gro
	
	sleep 10
	cp npt.top ../afe/gromacs.top
	cd ../
	echo "NPT equilibration complete."

	# AFE
	echo "Starting AFE calculation..."
	cd afe
	gmx grompp -f gromacs.mdp -c gromacs.gro -p gromacs.top -t ../npt/npt.cpt -o gromacs.tpr -maxwarn 1
	gmx mdrun -v -deffnm gromacs -gpu_id $gpu 

#	echo "$PWD"
	cd 	

done < "$unbound_data"


