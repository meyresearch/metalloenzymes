#!/bin/bash
#SBATCH -o PATH_TO_LOGS/%x_%a.slurm.out
#SBATCH -e PATH_TO_LOGS/%x_%a.slurm.err
#SBATCH -n N_TASKS
#SBATCH --gres=gpu:N_GPUS
#SBATCH --cpus-per-gpu N_CPUS
#SBATCH --mem MEMORY

start=`date +%s`

# get ligand pair and lambda string from input variables
transformation=$1
lambdastring="$2"

# set variables in meze.py
engine=ENGINE
repeats=N_REPEATS 
outputs_dir=OUTPUTS_DIR

# read in lambda string
INPUT_FILE_STREAM=',' read -a lambdas <<< "$lambdastring"

# use slurm array task id as lambda window index
id=$SLURM_ARRAY_TASK_ID
lambda=${lambdas[$id]}

for stage in "unbound" "bound"
do
	for (( i=1; i<=$repeats; i++)) 
	do
		minimisation_directory=$outputs_dir/${engine}_${i}/$transformation/$stage/minimisation/lambda_$lambda
		echo "Minimisation directory is: $minimisation_directory"
		lambda_directory=$outputs_dir/${engine}_${i}/$transformation/$stage/lambda_$lambda
		echo "Lambda directory is: $lambda_directory"
		
		echo "Using ${engine}_${i} for $transformation at $stage lambda $lambda"
		
		echo "Minimising..."
		cd $minimisation_directory
		echo "cd $minimisation_directory"
		if [[ $engine == *"SOMD"* ]]; then	
			# Minimise
			$BSS_HOME/somd-freenrg -C ./somd.cfg -l $lambda -c ./somd.rst7 -t ./somd.prm7 -m ./somd.pert -p CUDA 1> somd.out 2> somd.err
		
			# Copy files to lambda directory to be used for the AFE run
			cp $minimisation_directory/sim_restart.s3 $lambda_directory/sim_restart.s3
			cp $minimisation_directory/sim_restart.s3.previous $lambda_directory/sim_restart.s3.previous
			cp $minimisation_directory/SYSTEM.s3 $lambda_directory/SYSTEM.s3
			
			echo "Running AFE transformation..."
			cd $lambda_directory
		
			$BSS_HOME/somd-freenrg -C ./somd.cfg -l $lambda -c ./somd.rst7 -t ./somd.prm7 -m ./somd.pert -p CUDA 1> somd.out 2> somd.err
		else
			echo "Engine $engine is not supported yet."
		fi
	done
done

end=`date +%s`
runtime=$((end - start))

echo "Finished in $runtime seconds, or $((runtime/60)) minutes"
