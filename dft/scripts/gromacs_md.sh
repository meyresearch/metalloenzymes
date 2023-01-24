#!/bin/bash

ligands=( 1 3 4 6 7 9 10 12 13 16 )
ligands=( 1 )
cd ..
for i in ${ligands[@]}
do
	echo "ligand $i"
	cd ligand_$i/
	
	echo "Starting minimisation"
	mkdir 01_minimisation
	cp vim2_solv.gro 01_minimisation/
	cp vim2_solv.top 01_minimisation/
	cp ../scripts/min.mdp 01_minimisation/
	cd 01_minimisation/
	gmx grompp -f min.mdp -c vim2_solv.gro -p vim2_solv.top -o em.tpr > em.out
	gmx mdrun -v -deffnm em -nb gpu -gpu_id 0 
	sleep 10
	echo "10 0" | gmx energy -f em.edr -o potential.xvg 
	cd ..
	
	echo "Minimisation done. Starting 5 ns NVT equilibration."
	mkdir 02_nvt
	cd 02_nvt
	echo "1 q" | gmx make_ndx -f ../01_minimisation/em.gro -o posre.ndx
	echo "1" | gmx genrestr -f ../01_minimisation/vim2_solv.gro -n posre.ndx -o posre.itp
	python $HOME/projects/metalloenzymes/dft/scripts/edit_topology.py $i 
	cp ../../scripts/nvt.mdp .
	gmx grompp -f nvt.mdp -c ../01_minimisation/em.gro -r ../01_minimisation/em.gro -p ../01_minimisation/vim2_solv.top -o nvt.tpr > nvt.out	
	gmx mdrun -v -deffnm nvt -nb gpu -gpu_id 0
	sleep 10
	echo "16 0" | gmx energy -f nvt.edr -o temperature.xvg
	cd ..
	
	echo "NVT equilibration done. Starting 5 ns NPT equilibration."
	mkdir 03_npt
	cd 03_npt
	cp ../../scripts/npt.mdp .
        gmx grompp -f npt.mdp -c ../02_nvt/nvt.gro -r ../02_nvt/nvt.gro -t ../02_nvt/nvt.cpt -p ../01_minimisation/vim2_solv.top -o npt.tpr	
	gmx mdrun -v -deffnm npt -nb -gpu_id 0
	sleep 10
	echo "18 0" | gmx energy -f npt.edr -o pressure.xvg
	echo "24 0" | gmx energy -f npt.edr -o density.xvg
	echo "NPT equilibration done. Starting 50 ns MD production run."
	cd ..

	mkir 04_md
	cd 04_md
	cp ../../scripts/md.mdp .
	gmx grompp -f md.mdp -c ../03_npt/npt.gro -t ../03_npt/npt.cpt -p ../01_minimisation/vim2_solv.top -o md.tpr
	gmx mdrun -v -deffnm md -nb -gpu_id 0
	sleep 10
	cd ..
	echo "MD production run done for ligand $i."
	cd ..
done
