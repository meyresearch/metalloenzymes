#!/bin/bash
# -o controls output. %A is replaced with job ID and %a with array index
#SBATCH -o /home/jguven/projects/metalloenzymes/slurm_logs/somd1_afe.%A.%a.slurm.out
#SBATCH -e /home/jguven/projects/metalloenzymes//slurm_logs/somd1_afe.%A.%a.slurm.err
# Only launch a max of 1 task
#SBATCH -n 1
# allocate 1 gpu per job
#SBATCH --gres=gpu:1
# allocate 8 processors per gpu
#SBATCH --cpus-per-gpu 10
#SBATCH --mem 4096
# set job name
#SBATCH --job-name=somd_afe

export BSS_HOME="$HOME/Software/miniconda3/envs/bss-d/"
export OPENMM_PLUGIN_DIR="$BSS_HOME/lib/plugins/"
export OPEN_NUM_THREADS=1

date
start=`date +%s`

idx=$SLURM_ARRAY_TASK_ID

ligand_1=$1
ligand_2=$2
lig_1=`echo $ligand_1 | sed "s/ligand_/lig_/g"` 
lig_2=`echo $ligand_2 | sed "s/ligand_/lig_/g"`

engine=$3
lambdastring="$4"

IFS=' ' read -a lambdas <<< "$lambdastring"
lambda=${lambdas[$idx]}

log_dir=$HOME/projects/metalloenzymes/slurm_logs/
if [[ ! -d $log_dir ]]; then
	mkdir $log_dir
fi

for stage in "bound" "unbound"
do

    	min_dir=$HOME/projects/alchemistry/kpc2_deprotonated_ligand/outputs/${engine}_1/$lig_1~$lig_2/$stage/minimisation/lambda_$lambda
    	echo $min_dir
    	lambda_dir=$HOME/projects/alchemistry/kpc2_deprotonated_ligand/outputs/${engine}_1/$lig_1~$lig_2/$stage/lambda_$lambda
      
	echo "using $engine for $lig_1 and $lig_2, at lambda $lambda"
	echo "minimising"
      	echo $min_dir
	cd $min_dir
	if [[ $engine == *"SOMD"* ]]; then

		$BSS_HOME/bin/somd-freenrg -C ./somd.cfg -l $lambda -c ./somd.rst7 -t ./somd.prm7 -m ./somd.pert -p CUDA 1> somd.log 2> somd.err
		
		cp $min_dir/sim_restart.s3 $lambda_dir/sim_restart.s3
		cp $min_dir/sim_restart.s3.previous $lambda_dir/sim_restart.s3.previous
		cp $min_dir/SYSTEM.s3 $lambda_dir/SYSTEM.s3

		echo "perturbing"
		cd $lambda_dir	

		$BSS_HOME/bin/somd-freenrg -C ./somd.cfg -l $lambda -c ./somd.rst7 -t ./somd.prm7 -m ./somd.pert -p CUDA 1> somd.log 2> somd.err

	fi
		
done

end=`date +%s`
runtime=$((end - start))
echo "Finished in $runtime seconds, OR $((runtime/60)) minutes."

