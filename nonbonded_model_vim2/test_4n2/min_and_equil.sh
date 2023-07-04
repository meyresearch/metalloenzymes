#!/bin/bash
#SBATCH -o /home/jguven/projects/metalloenzymes/nonbonded_model_vim2/step_4n2/harmonic_restraints_md/slurm_log/min_and_equil.%A.%a.slurm.out
#SBATCH -e /home/jguven/projects/metalloenzymes/nonbonded_model_vim2/step_4n2/harmonic_restraints_md/slurm_log/min_and_equil.%A.%a.slurm.err
# Only launch a max of 1 task
#SBATCH -n 1
# allocate 1 gpu per job
#SBATCH --gres=gpu:1
# allocate 8 processors per gpu
#SBATCH --cpus-per-gpu 10
#SBATCH --mem 4096
# set job name
#SBATCH --job-name=min_and_equil


echo "starting 1min at $date"

$AMBERHOME/bin/pmemd.cuda -O -i 1min.in\
            -o 1min.out -p vim2_solv.prmtop -c vim2_solv.inpcrd -r 1min.rst7\
            -inf 1min.info -ref vim2_solv.inpcrd -x mdcrd.1min
echo "ending 1min at $date"
wait
echo "starting 2mdheat at $date"

$AMBERHOME/bin/pmemd.cuda -O -i 2mdheat.in\
            -o 2mdheat.out -p vim2_solv.prmtop -c 1min.rst7 -r 2mdheat.rst7\
            -inf 2mdheat.info -ref 1min.rst7 -x mdcrd.2mdheat
echo "ending 2mdheat at $date"

wait
echo "starting 3md at $date"

$AMBERHOME/bin/pmemd.cuda -O -i 3md.in\
            -o 3md.out -p vim2_solv.prmtop -c 2mdheat.rst7 -r 3md.rst7\
            -inf 3md.info -ref 2mdheat.rst7 -x mdcrd.3md
echo "ending 3md at $date"

wait
echo "starting 4md at $date"

$AMBERHOME/bin/pmemd.cuda -O -i 4md.in\
            -o 4md.out -p vim2_solv.prmtop -c 3md.rst7 -r 4md.rst7\
            -inf 4md.info -ref 3md.rst7 -x mdcrd.4md
echo "ending 4md at $date"

wait
echo "starting 5min at $date"

$AMBERHOME/bin/pmemd.cuda -O -i 5min.in\
            -o 5min.out -p vim2_solv.prmtop -c 4md.rst7 -r 5min.rst7\
            -inf 5min.info -ref 4md.rst7 -x mdcrd.5min
echo "ending 5min at $date"

wait
echo "starting 6md at $date"

$AMBERHOME/bin/pmemd.cuda -O -i 6md.in\
            -o 6md.out -p vim2_solv.prmtop -c 5min.rst7 -r 6md.rst7\
            -inf 6md.info -ref 5min.rst7 -x mdcrd.6md
echo "ending 6md at $date"

wait
echo "starting 7md at $date"

$AMBERHOME/bin/pmemd.cuda -O -i 7md.in\
            -o 7md.out -p vim2_solv.prmtop -c 6md.rst7 -r 7md.rst7\
            -inf 7md.info -ref 6md.rst7 -x mdcrd.7md
echo "ending 7md at $date"

wait
echo "starting 8md at $date"

$AMBERHOME/bin/pmemd.cuda -O -i 8md.in\
            -o 8md.out -p vim2_solv.prmtop -c 7md.rst7 -r 8md.rst7\
            -inf 8md.info -ref 7md.rst7 -x mdcrd.8md
echo "ending 8md at $date"

wait
echo "starting 9md at $date"

$AMBERHOME/bin/pmemd.cuda -O -i 9md.in\
            -o 9md.out -p vim2_solv.prmtop -c 8md.rst7 -r 9md.rst7\
            -inf 9md.info -ref 8md.rst7 -x mdcrd.9md
echo "ending 9md at $date"
