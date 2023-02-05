#!/bin/bash
# -o controls output. %A is replaced with job ID and %a with array index
#SBATCH -o /home/jguven/projects/metalloenzymes/dft/slurm_logs/md.%a.slurm.out
#SBATCH -e /home/jguven/projects/metalloenzymes/dft//slurm_logs/md.%a.slurm.err
# Only launch a max of 1 task
#SBATCH -n 1
# allocate 1 gpu per job
#SBATCH --gres=gpu:1
# allocate 8 processors per gpu
#SBATCH --cpus-per-gpu 10
#SBATCH --mem 4096
# set job name
#SBATCH --job-name=dft_md
# specify number of threads
export OMP_NUM_THREADS=8

# Set the ligand number as the array task id
ligand=$SLURM_ARRAY_TASK_ID

dft_dir=/home/jguven/projects/metalloenzymes/dft/

echo "-----------------------------------------------------------------"
echo "MINIMISATION $ligand"

# If minimisation folder exists, delete the old one
min_dir=$dft_dir/ligand_$ligand/01_minimisation/
if [[ -d $min_dir ]]; then
	rm -r $min_dir
fi
mkdir $min_dir

cd $min_dir
cp $dft_dir/ligand_$ligand/vim2_solv.gro $min_dir
cp $dft_dir/ligand_$ligand/vim2_solv.top $min_dir
cp $dft_dir/scripts/min.mdp $min_dir

gmx grompp -f min.mdp -c vim2_solv.gro -p vim2_solv.top -o em.tpr -maxwarn 1 
srun gmx mdrun -ntmpi 1 -v -deffnm em

cd ../
echo "-----------------------------------------------------------------"
echo "-----------------------------------------------------------------"
echo "NVT $ligand"
nvt_dir=$dft_dir/ligand_$ligand/02_nvt/

# If nvt folder exists, remove it 
if [[ -d $nvt_dir ]]; then
	rm -r $nvt_dir
fi
mkdir $nvt_dir
cd $nvt_dir
# Create position restraints index file
echo "1 q" | gmx make_ndx -f $min_dir/em.gro -o posre.ndx
# Create topology file for position restraints
echo "1" | gmx genrestr -f $min_dir/vim2_solv.gro -n posre.ndx -o posre.itp
# Edit topology to include position restraints
python $dft_dir/scripts/edit_topology.py $i 

cp $dft_dir/scripts/nvt.mdp .

gmx grompp -f nvt.mdp -c $min_dir/em.gro -r $min_dir/em.gro -p $min_dir/vim2_solv.top -o nvt.tpr -maxwarn 1
srun gmx mdrun -ntmpi 1 -v -deffnm nvt 
cd ../


echo "-----------------------------------------------------------------"
echo "-----------------------------------------------------------------"
echo "NPT $ligand"
npt_dir=$dft_dir/ligand_$ligand/03_npt/

# If npt folder exists, remove it 
if [[ -d $npt_dir ]]; then
	rm -r $npt_dir
fi
mkdir $npt_dir
cd $npt_dir
cp $dft_dir/scripts/npt.mdp .

gmx grompp -f npt.mdp -c $nvt_dir/nvt.gro -r $nvt_dir/nvt.gro -t $nvt_dir/nvt.cpt -p $min_dir/vim2_solv.top -o npt.tpr -maxwarn 1
srun gmx mdrun -ntmpi 1 -v -deffnm npt 
cd ../
echo "-----------------------------------------------------------------"
echo "-----------------------------------------------------------------"
echo "MD PRODUCTION $ligand"
md_dir=$dft_dir/ligand_$ligand/04_md/
if [[ -d $md_dir ]]; then
	rm -r $md_dir
fi
mkdir $md_dir
cd $md_dir
cp $dft_dir/scripts/md.mdp .

gmx grompp -f md.mdp -c $npt_dir/npt.gro -t $npt_dir/npt.cpt -p $min_dir/vim2_solv.top -o md.tpr
srun gmx mdrun -ntmpi 1 -v -deffnm md
cd ../
echo "-----------------------------------------------------------------"
