#! /bin/bash
#SBATCH -o /home/jguven/projects/metalloenzymes/nonbonded_model_vim2/step_4n2/harmonic_restraints_md/slurm_log/production_4n2.%A.%a.slurm.out
#SBATCH -e /home/jguven/projects/metalloenzymes/nonbonded_model_vim2/step_4n2/harmonic_restraints_md/slurm_log/production_4n2.%A.%a.slurm.err
# Only launch a max of 1 task
#SBATCH -n 1
# allocate 1 gpu per job
#SBATCH --gres=gpu:1
# allocate 8 processors per gpu
#SBATCH --cpus-per-gpu 10
#SBATCH --mem 4096
# set job name
#SBATCH --job-name=4n2

input_dir=$HOME/projects/metalloenzymes/nonbonded_model_vim2/step_4n2/harmonic_restraints_md/input_files/

$AMBERHOME/bin/pmemd.cuda -O -i md.in -p $input_dir/vim2_solv.prmtop -c 9md.rst7 -ref 9md.rst7 -o md.mdout -r vim2_md.rst7 -x vim2_md.nc
