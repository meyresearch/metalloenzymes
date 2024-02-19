# MetalloEnZymE parameterisation program (MEZE)

**Authors**: J. Jasmin Güven

# Table of Contents

[1. Prerequisites](#1-prerequisites)

[2. System setup](#2-system-setup)

[3. Prepare meze](#3-prepare-meze)

[4. Run meze AFE](#4-run-meze-afe)


# 1. Prerequisites 

Software requirements: 
- [BioSimSpace](https://biosimspace.openbiosim.org/install.html#easy-installation-run-in-a-conda-environment)

Useful software:
- [SLURM](https://slurm.schedmd.com/quickstart_admin.html)

# 2. System setup

- Start with a pdb file of the protein. You can save it in e.g. `inputs/protein/` in your project directory
- If your protein contains metal ions, save these in pdb format in the above directory
- Create ligand files in `*.sdf` or `*.mol2` formats and save in the same directory, e.g. `inputs/ligands/`
- If you want to keep crystallographic waters, save these in the protein directory in pdb format

# 3. Prepare meze

The workflow starts by setting up the project directory tree with the use of a python script called `prepare.py`. This script can be ran directly from the command line, or with the use of a bash script as:
```
$MEZEHOME/01_prepare.sh <path-to-input-files>/protein.pdb <ligand_charge> kpc2
```

For example, with KPC-2 with ligands whose net charge is -1:

```
$MEZEHOME/01_prepare.sh inputs/protein/kpc2.input.pdb -1 
```

Example output:
```
WARNING:root:Warning: importing 'simtk.openmm' is deprecated.  Import 'openmm' instead.
INFO:numexpr.utils:Note: NumExpr detected 32 cores but "NUMEXPR_MAX_THREADS" not set, so enforcing safe limit of 8.
INFO:numexpr.utils:NumExpr defaulting to 8 threads.

No separate water file detected in working directory. Continuing without it.

Successfully prepared meze network for AFE calculations.
Run scripts saved in: /home/jguven/projects/alchemistry/kpc2/partially_protonated_ligand///afe/
```

For help with `prepare.py` run: `python prepare.py -h`

The project directory tree will then look something like: 

```
/path/to/project/
    |--- inputs/                Input files (set-up by user)
    |   |--- ligands/           Ligand .sdf or .mol2 files (set-up by user)
    |   |--- protein/           Protein .pdb file (set-up by user)
    |--- outputs/               Folder for saving the results of AFE runs
    |--- afe/                   Input files for running the meze workflow and AFE runs
    |   |--- protocol.dat       Datafile with network information
    |   |--- ligands.dat        Datafile containing ligand names, e.g. ligand_1, ligand_2, ...
    |   |--- meze_network.csv   Dataframe of transformations from lomap, including ligand names, lomap scores and lambda windows
    |   |--- 02_add_water.sh    Bash script to carry out solvation of unbound and bound ligands with solvate.py
    |   |--- 03_heat_meze.sh    Bash script to carry out minimisation and equilibration of unbound and bound ligands with equilibrate.py
    |   |--- 04_meze.sh         Bash script to carry out AFE preparation with meze.py
    |   |--- run_SOMD.sh        Bash script to run an individual AFE transformation with a given engine, e.g. SOMD
    |   |--- slurm_all.sh       SLURM script which executes all of the above bash scripts
    |--- equilibration/         Folder for saving the minimisation and equilibration outputs for unbound and bound ligands
    |--- logs/                  Folder for keeping SLURM log files
```

# 4. Run meze AFE

To run the entire workflow using SLURM, simply run:

```
/path/to/project/afe/./slurm_all.sh
```

This executes all of the above scripts (02-04) using SLURM.

Example output:

```
dos2unix: converting file /home/jguven/projects/alchemistry/kpc2/partially_protonated_ligand///afe/ligands.dat to Unix format...
dos2unix: converting file /home/jguven/projects/alchemistry/kpc2/partially_protonated_ligand///afe//meze_network.csv to Unix format...
Adding water with slurm job 41654
Heating meze with slurm job 41655
Preparing AFE with slurm job 
Submitted AFE slurm job 41657
Submitted AFE slurm job 41658
Submitted AFE slurm job 41659
Submitted AFE slurm job 41663
Submitted AFE slurm job 41664
Submitted AFE slurm job 41665
Submitted AFE slurm job 41666
Submitted AFE slurm job 41667
Submitted AFE slurm job 41668
Submitted AFE slurm job 41669
Submitted AFE slurm job 41670
Submitted AFE slurm job 41671
Submitted AFE slurm job 41672
Submitted AFE slurm job 41673
Submitted AFE slurm job 41674
Submitted AFE slurm job 41675
Submitted AFE slurm job 41676
Submitted AFE slurm job 41677
Submitted AFE slurm job 41678
Submitted AFE slurm job 41679
Submitted AFE slurm job 41680
Submitted AFE slurm job 41681
```

To check the status of submitted jobs, use `squeue`.

The above project directory three will then be updated as:

```
/path/to/project/
    |--- inputs/                            Input files (set-up by user)
    |   |--- ligands/                       Ligand .sdf or .mol2 files (set-up by user)
    |   |--- protein/                       Protein .pdb file (set-up by user)
    |--- outputs/                           Folder for saving the results of AFE runs
    |   |--- SOMD_1/                        First AFE run directory
    |   |   |--- ligand_1~ligand_2/         AFE transformation of ligand_1 to ligand_2
    |   |   |   |--- bound/                 Bound stage
    |   |   |   |   |--- lambda_0.0000/     First lambda window 0.0
    |   |   |   |   |--- ...                [Rest of the lambda windows until lambda 1.0]
    |   |   |--- ...                        [Other transformations]
    |   |--- SOMD_2/                        Second AFE run repeat (depends on given number of repeats)
    |   |--- SOMD_3/                        Third AFE run repeat (depends on given number of repeats)
    |--- afe/                               Input files for running the meze workflow and AFE runs [same as above]
    |--- equilibration/                     Folder for saving the minimisation and equilibration outputs for unbound and bound ligands
    |   |--- bound/                         
    |   |   |--- ligand_1
    |   |   |   |--- min/                   Minimisation
    |   |   |   |--- bb_r_nvt/              Backbone-restrained NVT equilibration
    |   |   |   |--- r_nvt/                 Restrained NVT equilibration
    |   |   |   |--- nvt/                   NVT equilibration
    |   |   |   |--- r_npt/                 Restrained NPT equilibration
    |   |   |   |--- npt/                   NPT equilibration
    |   |   |--- ...
    |   |--- unbound/
    |   |   |--- ligand_1
    |   |   |   |--- min/                   Minimisation
    |   |   |   |--- r_nvt/                 Restrained NVT equilibration
    |   |   |   |--- nvt/                   NVT equilibration
    |   |   |   |--- r_npt/                 Restrained NPT equilibration
    |   |   |   |--- npt/                   NPT equilibration
    |   |   |--- ...    
    |--- logs/                              Folder for keeping SLURM log files
    |   |--- add_water_1.slurm.out          Solvation SLURM output for ligand 1 etc
    |   |--- add_water_1.slurm.err          Solvation SLURM error output for ligand 1 etc
    |   |--- heat_meze_1.slurm.out          Minimisation and equilibration SLURM output for ligand 1 etc
    |   |--- heat_meze_1.slurm.err          Minimisation and equilibration SLURM error output for ligand 1 etc
    |   |--- meze_1.slurm.out               AFE preparation SLURM output for ligand 1 etc
    |   |--- meze_1.slurm.err               AFE preparation SLURM error output for ligand 1 etc
    |   |--- ligand_1_ligand_2_1.slurm.out  AFE calculation SLURM output for the first lambda window for transformation ligand_1~ligand_2
    |   |--- ligand_1_ligand_2_1.slurm.err  AFE calculation SLURM error output for the first lambda window for transformation ligand_1~ligand_2
    |   |--- ...
    
    ```