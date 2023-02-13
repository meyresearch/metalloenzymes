#!/bin/bash
# -o controls output. %A is replaced with job ID and %a with array index
#SBATCH -o /home/jguven/projects/metalloenzymes/slurm_logs/somd_afe.%A.%a.slurm.out
#SBATCH -e /home/jguven/projects/metalloenzymes//slurm_logs/somd_afe.%A.%a.slurm.err
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

transformations=( lig_1~lig_2 lig_1~lig_4 lig_10~lig_4 lig_10~lig_8 lig_11~lig_15 lig_11~lig_5 lig_12~lig_14 lig_12~lig_15 lig_13~lig_6 lig_13~lig_7 lig_14~lig_15 lig_15~lig_5 lig_16~lig_9 lig_2~lig_4 lig_4~lig_6 lig_6~lig_7 lig_6~lig_8  )

transformation=${transformations[$idx]}

engine=$1
#lambdastring=$2

#IFS=',' read -r -a lambdas <<< "$lambdastring"
lambda=0.9000

log_dir=$HOME/projects/metalloenzymes/slurm_logs/
if [[ ! -d $log_dir ]]; then
	mkdir $log_dir
fi

for stage in "bound" "unbound"
do
	lambda_dir=$HOME/projects/metalloenzymes/kpc2/outputs/$engine/$transformation/$stage/lambda_$lambda
	cd $lambda_dir	
	echo "using $engine for $transformation, at lambda $lambda"
	
	if [[ $engine == *"SOMD"* ]]; then

		$BSS_HOME/bin/somd-freenrg -C ./somd.cfg -l $lambda -c ./somd.rst7 -t ./somd.prm7 -m ./somd.pert -p CUDA 1> somd.log 2> somd.err

	# add gromacs afe here

	fi
		
	
done

end=`date +%s`
runtime=$((end - start))
echo "Finished in $runtime seconds, OR $((runtime/60)) minutes."

