#!/bin/bash

INPUT_FILE=$1
LIGAND_CHARGE=$2

python $MEZEHOME/prepare.py --input-pdb-file $INPUT_FILE --ligand-charge $LIGAND_CHARGE
