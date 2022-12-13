
gpu=$1
bound_data=$2

while IFS= read -r bound_transformation
do
	cd $bound_transformation
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
	gmx grompp -f r_nvt.mdp -c r_nvt.gro -r r_nvt.gro -p r_nvt.top -o r_nvt.tpr
	gmx mdrun -v -deffnm r_nvt -gpu_id $gpu -c ../bb_r_nvt/bb_r_nvt.gro
	
	sleep 10
	cp r_nvt.top ../bb_r_nvt/bb_r_nvt.top
	cd ../
	echo "Restrained NVT complete."

	# Backbone-restrained NVT
	echo "Starting backbone-restrained NVT equilibration..."
	cd bb_r_nvt
	gmx grompp -f bb_r_nvt.mdp -c bb_r_nvt.gro -r bb_r_nvt.gro -t ../r_nvt/r_nvt.cpt -p bb_r_nvt.top -o bb_r_nvt.tpr
	gmx mdrun -v -deffnm bb_r_nvt -gpu_id $gpu -c ../nvt/nvt.gro

	sleep 10 
	cp bb_r_nvt.top ../nvt/nvt.top
	cd ../
	echo "Backbone-restrained NVT complete"

	# NVT
	echo "Starting NVT equilibration..."
	cd nvt
	gmx grompp -f nvt.mdp -c nvt.gro -t ../bb_r_nvt/bb_r_nvt.cpt -p nvt.top -o nvt.tpr
	gmx mdrun -v -deffnm nvt -gpu_id $gpu -c ../npt/npt.gro

	sleep 10
	cp nvt.top ../npt/npt.top
	cd ../
	echo "NVT equilibration complete."

	# NPT
	echo "Starting NPT equilibration..."
	cd npt
	gmx grompp -f npt.mdp -c npt.gro -p npt.top -t ../nvt/nvt.cpt -o npt.tpr
	gmx mdrun -v -deffnm npt -gpu_id $gpu -c ../afe/gromacs.gro
	
	sleep 10
	cp npt.top ../afe/gromacs.top
	cd ../
	echo "NPT equilibration complete."

	# AFE
	echo "Starting AFE calculation..."
	cd afe
	gmx grompp -f gromacs.mdp -c gromacs.gro -p gromacs.top -t ../npt/npt.cpt -o gromacs.tpr
	gmx mdrun -v -deffnm gromacs -gpu_id $gpu 

#	echo "$PWD"
	cd 	

done < "$bound_data"
