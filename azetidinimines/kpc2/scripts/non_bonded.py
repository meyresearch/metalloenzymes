import functions as fn
import pathlib
import glob
import os
import shutil


home_directory = str(pathlib.Path.home())
dft_directory = home_directory + "/projects/metalloenzymes/dft/"
nonbonded =  dft_directory + "nonbonded/" 

ligand_folders = sorted(glob.glob(dft_directory + "ligand_*"))
ligands = [path.split("/")[-1] for path in ligand_folders]
nonbonded_ligand_list = [nonbonded + ligand for ligand in ligands]

nonbonded_ligand_directories = [fn.create_dirs(path) + "/" for path in nonbonded_ligand_list]
amber_directories = [fn.create_dirs(path + "/amber/") for path in nonbonded_ligand_list]
min_directories = [fn.create_dirs(path + "/01_minimisation/") for path in amber_directories]
heat_directories = [fn.create_dirs(path + "/02_heat/") for path in amber_directories]
relax_directories = [fn.create_dirs(path + "/03_relax/") for path in amber_directories]
lower_directories = [fn.create_dirs(path + "/04_lower/") for path in amber_directories]
bb_min_directories = [fn.create_dirs(path + "/05_bb_min/") for path in amber_directories]
bb_relax_directories = [fn.create_dirs(path + "/06_bb_relax/") for path in amber_directories]
reduce_directories = [fn.create_dirs(path + "/07_reduce/") for path in amber_directories]
continue_directories = [fn.create_dirs(path + "/08_continue/") for path in amber_directories]
free_directories = [fn.create_dirs(path + "/09_free/") for path in amber_directories]
md_directories = [fn.create_dirs(path + "/10_md/") for path in amber_directories]

old_parmed_files = [sorted(glob.glob(path + "/parmed_*.out")) for path in ligand_folders]

for i in range(len(old_parmed_files)):
    ligand_parmed_list = old_parmed_files[i]
    for j in range(len(ligand_parmed_list)):
        parmed_file = ligand_parmed_list[j]
        filename = parmed_file.split("/")[-1]
        old_file = nonbonded_ligand_directories[i] + filename
        if not os.path.isfile(old_file):
            shutil.copyfile(parmed_file, old_file)

for i in range(len(nonbonded_ligand_directories)):
    data = []
    parmed_list = [nonbonded_ligand_directories[i] + "parmed_1.out",
                   nonbonded_ligand_directories[i] + "parmed_2.out"]

    for parmed in parmed_list:
        with open(parmed, "r") as file:
            lines = file.readlines()

        start_index = [i for i, string in enumerate(lines) if "Atom" in string][0] + 1
        end_index = start_index + 5

        for line in lines[start_index:end_index]:
            split_line = line.split()
            for string in split_line:
                try:
                    value = float(string)
                    data.append(value)
                except ValueError:
                    pass

    iat_1 = [int(data[i]) for i in range(0, len(data), 4)]
    iat_2 = [int(data[i]) for i in range(1, len(data), 4)]
    bond_length = [data[i] for i in range(2, len(data), 4)]
    force_constant = [data[i] for i in range(3, len(data), 4)]
    comment = "#Harmonic restraints for the restrained nonbonded model of VIM2\n"
    restraint_lines = [comment]

    for i in range(len(iat_1)):
        restraint = f"&rst iat={iat_1[i]},{iat_2[i]}, r1=0., r2={bond_length[i]}, r3={bond_length[i]}, r4=100., rk2={force_constant[i]}, rk3={force_constant[i]},/\n"
        restraint_lines.append(restraint)

    restraint_file = nonbonded_ligand_directories[i] + "disres.rst"
    with open(restraint_file, "w") as file:
        file.writelines(restraint_lines)
