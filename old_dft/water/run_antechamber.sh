#!/bin/bash

for i in $(seq 1 5)
do
	antechamber -fi pdb -fo mol2 -i WT$i.pdb -o WT$i.mol2 -at amber -c bcc -pf y
done
