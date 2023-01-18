#!/bin/bash

cd ../kpc2/outputs/GROMACS/

for transformation in */ 
do
	cd $transformation
	echo $transformation
        for stage in */
	do
		cd $stage
		for lambda in */
		do
			cd $lambda
			for dir in */
			do
				rm -r $dir/*			
			done
			cd ../
		done
		cd ../
	done
	cd ../	
done
