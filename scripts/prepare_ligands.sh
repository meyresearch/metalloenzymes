#!/bin/bash

ligands=( 1 2 3 4 5 6 7 9 10 11 12 13 14 15 16)


all_ligands="$HOME/projects/alchemistry/vim2_deprotonated_ligand/inputs/ligands"
parent_dir="$HOME/projects/alchemistry/vim2_deprotonated_ligand/parameterisation/"
# steps=( 1 2 )
# common_files="$parent_dir/common_files/*"
for i in ${ligands[@]}
do
	ligand_directory="$parent_dir/ligand_$i"
	if [[ ! -d $ligand_directory ]]; then
		mkdir $ligand_directory
	fi
	cd $ligand_directory
	# cp -r $common_files .
	cp $all_ligands/ligand_$i.pdb .
	echo $ligand_directory
	# for j in ${steps[@]}
	# do
		# step_directory="$ligand_directory/step_4n$j"
		# cd $step_directory
		# echo $step_directory
	grep -v '^CONECT' ligand_$i.pdb > L8J.pdb
	sed -i "s/HETATM/ATOM  /g" L8J.pdb
	sed -i "s/\*/\ /g" L8J.pdb
	antechamber -fi pdb -fo mol2 -i L8J.pdb -o L8J.mol2 -at gaff2 -c bcc -pf y -nc -2
	parmchk2 -i L8J.mol2 -o L8J.frcmod -f mol2
	sed -i -e :a -e '/^\n*$/{$d;N;ba' -e '}' L8J.frcmod
	sed -i -e '$a\  ' L8J.frcmod		
		# tail LIG.frcmod		
		# sed -i "s/L8J/LIG/g" vim2.in
	# done
done

