#!/bin/bash

sbatch --parsable --array=1-15 md_only.sh

