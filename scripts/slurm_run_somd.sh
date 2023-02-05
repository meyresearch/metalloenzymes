#!/bin/bash
#SBATCH -n 1 
#SBATCH --gres=gpu:0
#SBATCH --cpus-per-task=16
#SBATCH --job-name=somd-afe
#SBATCH -o ../slurm_logs/somd_afe_%A_%a.out
#SBATCH -e ../slurm_logs/somd_afe_%A_%a.err

# sys arg: $1: transformation, $2: number of windows


# Specify the number of lambda windows as the third sys arg
# TODO Should read these in from a file in future
if [ $2 = 11 ]; then
lambdas=( 0.0000 0.1000 0.2000 0.3000 0.4000 0.5000 0.6000 0.7000 0.8000 0.9000 1.0000 )
fi
if [ $2 = 17 ]; then
lambdas=( 0.0000 0.0625 0.1250 0.1875 0.2500 0.3125 0.3750 0.4375 0.5000 0.5625 0.6250 0.6875 0.7500 0.8125 0.8750 0.9375 1.0000 )
fi

lambda=${lambdas[SLURM_ARRAY_TASK_ID]}

cd ../outputs/SOMD/$1/
if [[ ! -d ../outputs/SOMD/$1/ ]]; then
	echo '../outputs/SOMD/$1 does not exist. AFE run aborted.'
	exit
fi

for directory in 'bound' 'free'
do
	for lambda in '${lambdas[@]}'
	do
		cd lambda_$lambda
#		somd-freenrg -C ./somd.cfg -l $lambda -d 0 -c ./somd.rst7 -t ./somd.prm7 -m ./somd.pert -p CUDA
		echo 'somd-freenrg -C ./somd.cfg -l $lambda -d 0 -c ./somd.rst7 -t ./somd.prm7 -m ./somd.pert -p CUDA'
		cd ../
	done
done

