#!/bin/bash

sbatch --parsable --array=1-15 run_plot.sh

