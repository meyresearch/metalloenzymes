#!/bin/bash

eval "$(conda shell.bash hook)"
conda activate obss


INPUT_FILE=$1
LIGAND_CHARGE=$2
#RUNTIME=$3

python $MEZEHOME/prepare.py --input-pdb-file $INPUT_FILE --ligand-charge $LIGAND_CHARGE 
