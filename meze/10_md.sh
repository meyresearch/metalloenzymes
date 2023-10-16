#!/bin/bash
#SBATCH -n 1
#SBATCH --gres=gpu:3
#SBATCH --cpus-per-gpu 10
#SBATCH --mem 4069

python $MEZEHOME/prepare_md_0.py "ligand_1" -if "$HOME/projects/alchemistry/vim2/deprotonated_ligand/md/inputs/protein/vim2.input.pdb" -c -2 -g vim2