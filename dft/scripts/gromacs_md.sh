#!/bin/bash
# -o controls output. %A is replaced with job ID and %a with array index
#SBATCH -o /home/jguven/projects/metalloenzymes/dft/slurm_logs/md.%A.%a.slurm.out
#SBATCH -e /home/jguven/projects/metalloenzymes/dft//slurm_logs/md.%A.%a.slurm.err
# Only launch a max of 1 task
#SBATCH -n 1
# allocate 1 gpu per job
#SBATCH --gres=gpu:1
# allocate 8 processors per gpu
#SBATCH --cpus-per-gpu 8
# set job name
#SBATCH --job-name=dft_md

# specify number of threads
export OMP_NUM_THREADS=8


ligands=( 1 3 4 6 7 9 10 12 13 16 )
ligands=( 1 )

# Create folder for slurm logs
if [[ ! -d /home/jguven/projects/metalloenzymes/dft//slurm_logs ]]; then
	mkdir /home/jguven/projects/metalloenzymes/dft/slurm_logs
fi

cd ..
for i in ${ligands[@]}
do
	echo "ligand $i"
	cd ligand_$i/
	
	echo "Starting minimisation"
	# If minimisation folder exists, remove it to get rid of old data (especially with testing)
	if [[ -d 01_minimisation ]]; then
		rm -r 01_minimisation
	fi
	# If minimisation folder doesn't exist, create one
	if [[ ! -d 01_minimisation ]]; then
		mkdir 01_minimisation
	fi
	# Copy gro and top files to minimisation folder
	cp vim2_solv.gro 01_minimisation/
	cp vim2_solv.top 01_minimisation/
	# Copy mdp file to minimisation folder
	cp ../scripts/min.mdp 01_minimisation/
	cd 01_minimisation/
	# grompp
	gmx grompp -f min.mdp -c vim2_solv.gro -p vim2_solv.top -o em.tpr > em.out
	# mdrun
	srun srun gmx mdrun -ntmpi 1 -v -deffnm em  
	# wait
	sleep 10
	# get potential energy during minimisation
	#echo "10 0" | gmx energy -f em.edr -o potential.xvg 
	cd ..
	
	echo "Minimisation done. Starting 5 ns NVT equilibration."
	# If nvt folder exists, remove it 
	if [[ -d 02_nvt ]]; then
		rm -r 02_nvt
	fi
	# If nvt folder doesn't exist, create one
	if [[ ! -d 02_nvt ]]; then
		mkdir 02_nvt
	fi
	cd 02_nvt
	# Create position restraints index file
	#echo "1 q" | gmx make_ndx -f ../01_minimisation/em.gro -o posre.ndx
	# Create topology file for position restraints
	#echo "1" | gmx genrestr -f ../01_minimisation/vim2_solv.gro -n posre.ndx -o posre.itp
	# Edit topology to include position restraints
	python $HOME/projects/metalloenzymes/dft/scripts/edit_topology.py $i 
	# Copy mdp file to nvt folder
	cp ../../scripts/nvt.mdp .
	# grompp
	gmx grompp -f nvt.mdp -c ../01_minimisation/em.gro -r ../01_minimisation/em.gro -p ../01_minimisation/vim2_solv.top -o nvt.tpr > nvt.out	
	# mdrun
	srun gmx mdrun -ntmpi 1 -v -deffnm nvt 
	# wait
	sleep 10
	# get temperature during nvt
	#echo "16 0" | gmx energy -f nvt.edr -o temperature.xvg
	cd ..
	
	echo "NVT equilibration done. Starting 5 ns NPT equilibration."
	# if npt folder exists, remove it
	if [[ -d 03_npt ]]; then
		rm -r 03_npt
	fi
	# if npt folder doesn't exist, make one
	if [[ ! -d 03_npt ]]; then
		mkdir 03_npt
	fi
	cd 03_npt
	# copy mdp file to npt folder
	cp ../../scripts/npt.mdp .
	# grompp
        gmx grompp -f npt.mdp -c ../02_nvt/nvt.gro -r ../02_nvt/nvt.gro -t ../02_nvt/nvt.cpt -p ../01_minimisation/vim2_solv.top -o npt.tpr	
	# mdrun
	srun gmx mdrun -ntmpi 1 -v -deffnm npt 
	# wait
	sleep 10
	# get pressure and density during npt
	#echo "18 0" | gmx energy -f npt.edr -o pressure.xvg
	#echo "24 0" | gmx energy -f npt.edr -o density.xvg
	echo "NPT equilibration done. Starting 50 ns MD production run."
	cd ..
	# if md folder exists, remove it
	if [[ -d 04_md ]]; then
		rm -r 04_md
	fi
	# if md folder doesn't exist, create one
	if [[ ! -d 04_md ]]; then
		mkir 04_md
	fi
	cd 04_md
	# copy mdp file to md folder	
	cp ../../scripts/md.mdp .
	# grompp
	gmx grompp -f md.mdp -c ../03_npt/npt.gro -t ../03_npt/npt.cpt -p ../01_minimisation/vim2_solv.top -o md.tpr
	# mdrun
	srun gmx mdrun -ntmpi 1 -v -deffnm md 
	# wait
	sleep 10
	cd ..
	echo "MD production run done for ligand $i."
	cd ..
done
