#!/bin/bash

cd ../kpc2/outputs/SOMD/

for transformation in */ 
do
	cd $transformation
	echo $transformation
        for stage in */
	do
		cd $stage
		for lambda in */
		do
			rm -r $lambda/*			
		done
		cd ../
	done
	cd ../	
done
