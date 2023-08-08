#!/bin/bash

export AMBERHOME="/home/jguven/Software/amber22/"

all_ligands=$1
parent_dir=$2
charge=$3
ligand=$4
for i in {1..16}
do
	ligand_directory="$parent_dir/ligand_$i"
	if [[ ! -d $ligand_directory ]]; then
		mkdir $ligand_directory
	fi
	cd $ligand_directory
	cp $all_ligands/ligand_$i.pdb .
	echo $ligand_directory
	grep -v '^CONECT' ligand_$i.pdb > LIG.pdb
	sed -i "s/$ligand/LIG/g" LIG.pdb
	sed -i "s/HETATM/ATOM  /g" LIG.pdb
	sed -i "s/\*/\ /g" LIG.pdb
	antechamber -fi pdb -fo mol2 -i LIG.pdb -o LIG.mol2 -at gaff2 -c bcc -pf y -nc $charge
	parmchk2 -i LIG.mol2 -o LIG.frcmod -f mol2
	sed -i -e :a -e '/^\n*$/{$d;N;ba' -e '}' LIG.frcmod
	sed -i -e '$a\  ' LIG.frcmod		
done

