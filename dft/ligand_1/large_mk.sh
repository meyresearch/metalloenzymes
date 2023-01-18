#!/bin/bash
#$ -N mk_1
#$ -cwd
#$ -l h_rt=48:00:00
#$ -l h_vmem=2G
#$ -pe sharedmem 8
#  job name: -N
#  use the current working directory: -cwd
#  runtime limit of 24 hours: -l h_rt
#  memory limit of 2 Gbyte per slot: -l h_vmem
#  parallel environment and no. of slots: -pe

# Initialise the environment modules
. /etc/profile.d/modules.sh

# Export environment variables
export g16root=/exports/applications/apps/community/chem
export GAUSS_SCRDIR=$TMPDIR
source $g16root/g16/bsd/g16.profile

/exports/applications/apps/community/chem/g16/g16 vim2_large_mk.com
