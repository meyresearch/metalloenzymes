#!/bin/bash

cd ../kpc2/outputs/GROMACS/lig_2~lig_4/

for i in */ 
do
	cd unbound/
	for j in */
	do	
		rm -r $j	
	done
	cd ../
	cd bound/
	for j in */
	do	
		rm -r $j
	done
	cd ../
done
