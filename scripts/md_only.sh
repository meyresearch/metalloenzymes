#!/bin/bash
#SBATCH -o /home/jguven/projects/metalloenzymes/slurm_logs/step_4n2_min_and_equil.%A.%a.slurm.out
#SBATCH -e /home/jguven/projects/metalloenzymes/slurm_logs/step_4n2_min_and_equil.%A.%a.slurm.err
# Only launch a max of 1 task
#SBATCH -n 1
# allocate 1 gpu per job
#SBATCH --gres=gpu:1
# allocate 8 processors per gpu
#SBATCH --cpus-per-gpu 10
#SBATCH --mem 4096
# set job name
#SBATCH --job-name=min_and_equil
date

idx=$SLURM_ARRAY_TASK_ID
input_dir=$HOME/projects/metalloenzymes/nonbonded_model_vim2/ligand_$idx/
# nb_dir=$HOME/projects/metalloenzymes/nonbonded_model_vim2/step_4n2/harmonic_restraints_md/

cd $input_dir
echo $input_dir

production_dir=$input_dir/production/
# if [[ ! -d $production_dir ]]; then
#     mkdir $production_dir
# fi

cd $production_dir

# cp $nb_dir/*.in .
cp $HOME/projects/metalloenzymes/nonbonded_model_vim2/common_files/production/md.in .

echo "starting md at $date"

$AMBERHOME/bin/pmemd.cuda -O -i md.in -p vim2_solv.prmtop -c 9md.rst7 -ref 9md.rst7 -o md.mdout -r vim2_md.rst7 -x vim2_md.nc

echo "ending md at $date"
