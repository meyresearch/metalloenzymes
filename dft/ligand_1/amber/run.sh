#!/bin/bash
# -o controls output. %A is replaced with job ID and %a with array index
#SBATCH -o /home/jguven/projects/metalloenzymes/dft/slurm_logs/nonbonded_md.%A.%a.slurm.out
#SBATCH -e /home/jguven/projects/metalloenzymes/dft/slurm_logs/nonbonded_md.%A.%a.slurm.err
# Only launch a max of 1 task
#SBATCH -n 1
# allocate 1 gpu per job
#SBATCH --gres=gpu:1
# allocate 8 processors per gpu
#SBATCH --cpus-per-gpu 10
#SBATCH --mem 4096
# set job name
#SBATCH --job-name=nonbonded_md

run_dir=$HOME/projects/metalloenzymes/dft/ligand_1/amber/

min_dir=$run_dir/01_minimisation/
echo $min_dir
cd $min_dir
srun pmemd.cuda -O -i min.in -o min.out -p vim2_solv.prmtop -c vim2_solv.inpcrd -r min.rst7 -inf min.info -ref vim2_solv.inpcrd -x mdcrd.min

heat_dir=$run_dir/02_heat/
echo $heat_dir
cd $heat_dir
srun pmemd.cuda -O -i heat.in -o heat.out -p $min_dir/vim2_solv.prmtop -c $min_dir/min.rst7 -r heat.rst7 -inf heat.info -ref $min_dir/min.rst7 -x mdcrd.heat

relax_dir=$run_dir/03_relax/
echo $relax_dir
cd $relax_dir
srun pmemd.cuda -O -i relax.in -o relax.out -p $min_dir/vim2_solv.prmtop -c $heat_dir/heat.rst7 -r relax.rst7 -inf relax.info -ref $heat_dir/heat.rst7 -x mdcrd.relax

lower_dir=$run_dir/04_lower/
echo $lower_dir
cd $lower_dir
srun pmemd.cuda -O -i lower.in -o lower.out -p $min_dir/vim2_solv.prmtop -c $relax_dir/relax.rst7 -r lower.rst7 -inf lower.info -ref $relax_dir/relax.rst7 -x mdcrd.lower

bb_min_dir=$run_dir/05_bb_min/
echo $bb_min_dir
cd $bb_min_dir
srun pmemd.cuda -O -i bb_min.in -o bb_min.out -p $min_dir/vim2_solv.prmtop -c $lower_dir/lower.rst7 -r bb_min.rst7 -inf bb_min.info -ref $lower_dir/lower.rst7 -x mdcrd.bb_min

bb_relax_dir=$run_dir/06_bb_relax/
echo $bb_relax_dir
cd $bb_relax_dir
srun pmemd.cuda -O -i bb_relax.in -o bb_relax.out -p $min_dir/vim2_solv.prmtop -c $bb_min_dir/bb_min.rst7 -r bb_relax.rst7 -inf bb_relax.info -ref $bb_min_dir/bb_min.rst7 -x mdcrd.bb_relax

reduce_dir=$run_dir/07_reduce/
echo $reduce_dir
cd $reduce_dir
srun pmemd.cuda -O -i reduce.in -o reduce.out -p $min_dir/vim2_solv.prmtop -c $bb_relax_dir/bb_relax.rst7 -r reduce.rst7 -inf reduce.info -ref $bb_relax_dir/bb_relax.rst7 -x mdcrd.reduce

continue_dir=$run_dir/08_continue/
echo $continue_dir
cd $continue_dir
srun pmemd.cuda -O -i continue.in -o continue.out -p $min_dir/vim2_solv.prmtop -c $reduce_dir/reduce.rst7 -r continue.rst7 -inf continue.info -ref $reduce_dir/reduce.rst7 -x mdcrd.continue

free_dir=$run_dir/09_free/
echo $free_dir
cd $free_dir
srun pmemd.cuda -O -i free.in -o free.out -p $min_dir/vim2_solv.prmtop -c $continue_dir/continue.rst7 -r free.rst7 -inf free.info -ref $continue_dir/continue.rst7 -x mdcrd.free

md_dir=$run_dir/10_md/
cp $min_dir/vim2_solv.prmtop $md_dir/vim2_equil.prmtop
cp $free_dir/free.rst7 $md_dir/vim2_equil.rst7

cd $md_dir
echo $md_dir
srun pmemd.cuda -O -i md.in -o md.out -p vim2_equil.prmtop -c vim2_equil.rst7 -r md.rst7 -inf md.info -ref vim2_equil.rst7 -x md.nc
