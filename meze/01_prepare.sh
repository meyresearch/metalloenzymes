#!/bin/bash


INPUT_FILE=$1
LIGAND_CHARGE=$2
NAME=$3

#HAVE TO ADD --non-metal somehow?

python $MEZEHOME/prepare.py --input-pdb-file $INPUT_FILE --ligand-charge $LIGAND_CHARGE --group-name $NAME  
