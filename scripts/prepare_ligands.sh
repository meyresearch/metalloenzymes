#!/bin/bash

ligands=( 1 2 3 4 5 6 7 9 10 11 12 13 14 15 )

parent_dir="$HOME/projects/metalloenzymes/nonbonded_model_vim2/"
steps=( 1 2 )

for i in ${ligands[@]}
do
	ligand_directory="$parent_dir/ligand_$i"
	cd $ligand_directory
	echo $ligand_directory
	for j in ${steps[@]}
	do
		step_directory="$ligand_directory/step_4n$j"
		cd $step_directory
		echo $step_directory
		grep -v '^CONECT' ligand_$i.pdb > LIG.pdb
		sed -i "s/HETATM/ATOM  /g" LIG.pdb
		sed -i "s/L8J/LIG/g" LIG.pdb
		sed -i "s/MOL/LIG/g" LIG.pdb
		sed -i "s/\*/\ /g" LIG.pdb
		antechamber -fi pdb -fo mol2 -i LIG.pdb -o LIG.mol2 -at gaff2 -c bcc -pf y -nc -1
		# Convert re-numbered mol2 file to pdb 
		# obabel -i mol2 LIG.mol2 -o pdb -O LIG.pdb
		parmchk2 -i LIG.mol2 -o LIG.frcmod -f mol2
		sed -i -e :a -e '/^\n*$/{$d;N;ba' -e '}' LIG.frcmod
		sed -i -e '$a\  ' LIG.frcmod		
		# tail LIG.frcmod		
		# sed -i "s/L8J/LIG/g" vim2.in
	done
done

