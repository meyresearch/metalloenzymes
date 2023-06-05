#!/bin/bash
# -o controls output. %A is replaced with job ID and %a with array index
#SBATCH -o /home/jguven/projects/metalloenzymes/dft/slurm_logs/amber_md.%a.slurm.out
#SBATCH -e /home/jguven/projects/metalloenzymes/dft//slurm_logs/amber_md.%a.slurm.err
# Only launch a max of 1 task
#SBATCH -n 1
# allocate 1 gpu per job
#SBATCH --gres=gpu:1
# allocate 8 processors per gpu
#SBATCH --cpus-per-gpu 10
#SBATCH --mem 4096
# set job name
#SBATCH --job-name=dft_md

# Set the ligand number as the array task id
ligand=$SLURM_ARRAY_TASK_ID

dft_dir=/home/jguven/projects/metalloenzymes/dft/

echo "-----------------------------------------------------------------"
echo "MINIMISATION $ligand"

# If minimisation folder exists, delete the old one
min_dir=$dft_dir/ligand_$ligand/amber/01_minimisation/
if [[ -d $min_dir ]]; then
	rm -r $min_dir
fi
mkdir -p $min_dir

cd $min_dir
cp $dft_dir/ligand_$ligand/vim2_solv.inpcrd $min_dir
cp $dft_dir/ligand_$ligand/vim2_solv.prmtop $min_dir
cp $dft_dir/scripts/min.in $min_dir

srun pmemd.cuda -O -i min.in -o min.out -p vim2_solv.prmtop -c vim2_solv.inpcrd -r min.rst7 -inf min.info -ref vim2_solv.inpcrd -x mdcrd.min 

echo "-----------------------------------------------------------------"
echo "-----------------------------------------------------------------"
echo "HEAT $ligand"
heat_dir=$dft_dir/ligand_$ligand/amber/02_heat/

# If heat folder exists, remove it 
if [[ -d $heat_dir ]]; then
	rm -r $heat_dir
fi
mkdir -p $heat_dir
cd $heat_dir
cp $dft_dir/scripts/heat.in .

srun pmemd.cuda -O -i heat.in -o heat.out -p $min_dir/vim2_solv.prmtop -c $min_dir/min.rst7 -r heat.rst7 -inf heat.info -ref $min_dir/min.rst7 -x mdcrd.heat 




echo "-----------------------------------------------------------------"
echo "-----------------------------------------------------------------"
echo "RELAX $ligand"
relax_dir=$dft_dir/ligand_$ligand/amber/03_relax/

# If relax folder exists, remove it 
if [[ -d $relax_dir ]]; then
	rm -r $relax_dir
fi
mkdir -p $relax_dir
cd $relax_dir
cp $dft_dir/scripts/relax.in .

srun pmemd.cuda -O -i relax.in -o relax.out -p $min_dir/vim2_solv.prmtop -c $heat_dir/heat.rst7 -r relax.rst7 -inf relax.info -ref $heat_dir/heat.rst7 -x mdcrd.relax


echo "-----------------------------------------------------------------"
echo "-----------------------------------------------------------------"
echo "LOWER $ligand"
lower_dir=$dft_dir/ligand_$ligand/amber/04_lower/

# If lower folder exists, remove it 
if [[ -d $lower_dir ]]; then
	rm -r $lower_dir
fi
mkdir -p $lower_dir
cd $lower_dir
cp $dft_dir/scripts/lower.in .

srun pmemd.cuda -O -i lower.in -o lower.out -p $min_dir/vim2_solv.prmtop -c $relax_dir/relax.rst7 -r lower.rst7 -inf lower.info -ref $relax_dir/relax.rst7 -x mdcrd.lower


echo "-----------------------------------------------------------------"
echo "-----------------------------------------------------------------"
echo "BB MINIMISATION $ligand"
bb_min_dir=$dft_dir/ligand_$ligand/amber/05_bb_min/

# If bb min folder exists, remove it
if [[ -d $bb_min_dir ]]; then
	rm -r $bb_min_dir
fi
mkdir -p $bb_min_dir
cd $bb_min_dir
cp $dft_dir/scripts/bb_min.in .

srun pmemd.cuda -O -i bb_min.in -o bb_min.out -p $min_dir/vim2_solv.prmtop -c $lower_dir/lower.rst7 -r bb_min.rst7 -inf bb_min.info -ref $lower_dir/lower.rst7 -x mdcrd.bb_min


echo "-----------------------------------------------------------------"
echo "-----------------------------------------------------------------"
echo "BB RELAX $ligand"
bb_relax_dir=$dft_dir/ligand_$ligand/amber/06_bb_relax/

# If bb relax folder exists, remove it
if [[ -d $bb_relax_dir ]]; then
	rm -r $bb_relax_dir
fi
mkdir -p $bb_relax_dir
cd $bb_relax_dir
cp $dft_dir/scripts/bb_relax.in .

srun pmemd.cuda -O -i bb_relax.in -o bb_relax.out -p $min_dir/vim2_solv.prmtop -c $bb_min_dir/bb_min.rst7 -r bb_relax.rst7 -inf bb_relax.info -ref $bb_min_dir/bb_min.rst7 -x mdcrd.bb_relax


echo "-----------------------------------------------------------------"
echo "-----------------------------------------------------------------"
echo "REDUCE $ligand"
reduce_dir=$dft_dir/ligand_$ligand/amber/07_reduce/

# If reduce folder exists, remove it
if [[ -d $reduce_dir ]]; then
	rm -r $reduce_dir
fi
mkdir -p $reduce_dir
cd $reduce_dir
cp $dft_dir/scripts/reduce.in .

srun pmemd.cuda -O -i reduce.in -o reduce.out -p $min_dir/vim2_solv.prmtop -c $bb_relax_dir/bb_relax.rst7 -r reduce.rst7 -inf reduce.info -ref $bb_relax_dir/bb_relax.rst7 -x mdcrd.reduce


echo "-----------------------------------------------------------------"
echo "-----------------------------------------------------------------"
echo "CONTINUE $ligand"
continue_dir=$dft_dir/ligand_$ligand/amber/08_continue/

# If continue folder exists, remove it
if [[ -d $continue_dir ]]; then
	rm -r $continue_dir
fi
mkdir -p $continue_dir
cd $continue_dir
cp $dft_dir/scripts/continue.in .

srun pmemd.cuda -O -i continue.in -o continue.out -p $min_dir/vim2_solv.prmtop -c $reduce_dir/reduce.rst7 -r continue.rst7 -inf continue.info -ref $reduce_dir/reduce.rst7 -x mdcrd.continue


echo "-----------------------------------------------------------------"
echo "-----------------------------------------------------------------"
echo "FREE $ligand"
free_dir=$dft_dir/ligand_$ligand/amber/09_free/

# If free folder exists, remove it
if [[ -d $free_dir ]]; then
	rm -r $free_dir
fi
mkdir -p $free_dir
cd $free_dir
cp $dft_dir/scripts/free.in .

srun pmemd.cuda -O -i free.in -o free.out -p $min_dir/vim2_solv.prmtop -c $continue_dir/continue.rst7 -r free.rst7 -inf free.info -ref $continue_dir/continue.rst7 -x mdcrd.free


echo "-----------------------------------------------------------------"
echo "-----------------------------------------------------------------"
echo "MD PRODUCTION $ligand"
md_dir=$dft_dir/ligand_$ligand/amber/10_md/

if [[ -d $md_dir ]]; then
	rm -r $md_dir
fi

mkdir -p $md_dir
cd $md_dir
cp $dft_dir/scripts/md.mdp .
cp $min_dir/vim2_solv.prmtop vim2_equil.prmtop
cp $free_dir/free.rst7 vim2_equil.rst7
srun pmemd.cuda -O -i md.in -o md.out -p vim2_equil.prmtop -c vim2_equil.rst7 -inf md.info -ref vim2_equil.rst7 -r vim2_equil.rst7 -x md.nc


echo "-----------------------------------------------------------------"
