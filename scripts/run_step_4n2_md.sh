#!/bin/bash

sbatch --parsable --array=1-15 step_4n2_md.sh

