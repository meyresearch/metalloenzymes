
#!/bin/bash

ligands=( 1 2 3 4 5 6 7 9 10 11 12 13 14 15 16 )

cd ..
for i in ${ligands[@]}
do
	cd ligand_$i
	echo "ligand_$i"
	rsync -av pluto:/home/jguven/projects/metalloenzymes/dft/ligand_$i/vim2_solv.prmtop ../structures/lig_$i.prmtop
	rsync -av pluto:/home/jguven/projects/metalloenzymes/dft/ligand_$i/vim2_solv.inpcrd ../structures/lig_$i.inpcrd
	cd ../
done
