#!/bin/bash

ligands=( 1 2 3 4 5 6 7 9 10 11 12 13 14 15 16 )

cd ..
for i in ${ligands[@]}
do
	cd ligand_$i
	cp ligand_$i.pdb LIG_not_renamed_atoms.pdb
	sed -i "s/TWB/LIG/g" LIG_not_renamed_atoms.pdb
	sed -i "s/MOL/LIG/g" LIG_not_renamed_atoms.pdb
	sed -i "s/\*/\ /g" LIG_not_renamed_atoms.pdb
	antechamber -fi pdb -fo mol2 -i LIG_not_renamed_atoms.pdb -o LIG.mol2 -at gaff2 -c bcc -pf y -nc -1
	# Convert re-numbered mol2 file to pdb 
	obabel -i mol2 LIG.mol2 -o pdb -O LIG.pdb
	parmchk2 -i LIG.mol2 -o LIG.frcmod -f mol2
	#sed "$(( $(wc -l < LIG.frcmod) - 2 + 1 )), $ d" LIG.frcmod
	sed -i -e :a -e '/^\n*$/{$d;N;ba' -e '}' LIG.frcmod
	sed -i -e '$a\  ' LIG.frcmod		
	tail LIG.frcmod		
	#sed -i -e :a -e '/^\n*$/{$d;N;ba' -e '}' LIG.frcmod
	cd ..

done

