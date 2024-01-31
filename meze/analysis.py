import math
import shutil
import numpy as np
import argparse
import functions
import os
import subprocess as sp
import MDAnalysis as mda
import MDAnalysis.analysis.rms
import warnings
import logging
import meze
logger = logging.getLogger()
logger.setLevel(logging.CRITICAL)


def read_mbar(mbar_textfile):
    """
    Read in mbar.txt and return binding free energy and an error estimate from MBART

    Parameters:
    -----------
    mbar_textfile: 
        full path to mbar.txt

    Return:
    -------
    free_energy, error: tuple
        values for free energy and its associated error
    """
    with open(mbar_textfile, "r") as file:
        lines = file.readlines()
    result_line = 0
    for i, line in enumerate(lines):
        if "#MBAR free energy difference in kcal/mol" in line:
            result_line = i + 1
    free_energy, error = None, None
    try:
        if "#" in lines[result_line]:
            split_line = lines[result_line].split("#")
            warnings.warn(split_line[-1])
            try:
                free_energy = float(split_line[0].split(",")[0])
                error = float(split_line[0].split(",")[-1])
            except ValueError as error_message:
                print(f"Error: {error_message}")
    except IndexError as e:
        print(f"{e}: {mbar_textfile}")

    else:
        try:
            free_energy = float(lines[result_line].split(",")[0])
            error = float(lines[result_line].split(",")[-1])
        except ValueError as error_message:
            print(f"Error: {error_message}")  
   
    return free_energy, error 


def write_header(simfile, correct_header):
    """
    Take the new, correct header and append it to the start of the simfile

    Parameters:
    -----------
    simfile: str
        full path to the simfile in the specific lambda window
    correct_header: list
        list of the correct header lines for the specific lambda window

    Return:
    -------
    """
    data = []
    with open(simfile, "r") as file:

        for line in file:
            if "#" not in line and not line.isspace():
                data.append(line)

        # Check last value of datafile
        # If it's not the last frame (#TODO remove hard-coding this), delete the line
        # This is because there seem to be some random numbers at the end of the files and I'm not sure where they are coming from 
        
        for i, line in enumerate(data):
            split_line = line.split()
            if "#" not in line and len(split_line) < 16: 
                del data[i]

        file.seek(0, 0)
        header_and_data = correct_header + data
    with open(simfile, "w") as file:
        file.writelines(header_and_data)


def fix_simfile(protocol, transformation): 
    """
    Get the headers from minimisation simfile for each lambda file, which contain the correct headers required by MBAR.
    Nothing will be done if minimisation directories have already been moved.
    Note: Only for SOMD

    Parameters:
    -----------
    protocol: dict
        the protocol file as a dictionary
    transformation: str
        name of the transformation

    Return:
    -------
    """
    template = os.environ["MEZEHOME"] + "simfile_header.txt"
    with open(template, "r") as file:
        template_header = file.readlines()

    outputs = functions.path_exists(protocol["outputs"])
    engine = protocol["engine"]
    n_repeats = functions.check_int(protocol["repeats"])

    for i in range(1, n_repeats + 1): 

        path = f"{outputs}/{engine}_{i}/{transformation}/"

        unbound_directory = path + "unbound/"
        bound_directory = path + "bound/"

        unbound_minimisation_simfiles = functions.get_files(unbound_directory + "/minimisation/lambda_*/simfile.dat")
        unbound_simfiles = functions.get_files(unbound_directory + "lambda_*/simfile.dat")
        bound_simfiles = functions.get_files(bound_directory + "lambda_*/simfile.dat")
        
        old_unbound_files = [filename.replace(".dat", "_original") for filename in unbound_simfiles]
        old_bound_files = [filename.replace(".dat", "_original") for filename in bound_simfiles]

        _ = [shutil.copy(unbound_simfiles[i], old_unbound_files[i]) for i in range(len(old_unbound_files))]
        _ = [shutil.copy(bound_simfiles[i], old_bound_files[i]) for i in range(len(old_bound_files))]

        for i in range(len(unbound_minimisation_simfiles)):
            with open(unbound_minimisation_simfiles[i], "r") as file:
                minimisation_simfile = file.readlines()
            if "#" in minimisation_simfile[0]:
                with open(unbound_minimisation_simfiles[i], "r") as file:
                    for j, line in enumerate(file):
                        if "lambda" in line:
                            start = j
                        if "#" not in line:
                            end = j
                            break

                lambda_lines = minimisation_simfile[start:end]
                correct_lambda_header = template_header + lambda_lines

                write_header(unbound_simfiles[i], correct_lambda_header)
                write_header(bound_simfiles[i], correct_lambda_header)
            else: 
                continue
                

def save_results(protocol, transformation):
    """
    Using the protocol, get results from each repeat and write to file.

    Adapted from https://github.com/OpenBioSim/biosimspace_tutorials/blob/main/04_fep/02_RBFE/scripts/analysis.py 

    Parameters:
    -----------
    protocol: dict
        protocol file as a dictionary
    transformation: str
        name of the transformation

    Return:
    -------

    """
    engine = protocol["engine"]
    outputs = protocol["outputs"]
    n_repeats = functions.check_int(protocol["repeats"])

    for i in range(1, n_repeats + 1):
        path = f"{outputs}/{engine}_{i}/{transformation}/" 
        unbound = path + "unbound/"
        bound = path + "bound/"

        # Move minimisation directory out of transformation directory to avoid BioSimSpace error
        try:
            unbound_minimisation_directory = f"{outputs}/{engine}_{i}/" + "unbound_minimisations"
            bound_minimisation_directory =  f"{outputs}/{engine}_{i}/" + "bound_minimisations"
            functions.mkdir(unbound_minimisation_directory)
            functions.mkdir(bound_minimisation_directory)
            target_path = unbound_minimisation_directory + "/" + transformation + "/"
            source_path = unbound + "minimisation"
            shutil.move(source_path, target_path)
            target_path = bound_minimisation_directory + "/" + transformation + "/"
            source_path = bound + "minimisation"
            shutil.move(source_path, target_path)

        except FileNotFoundError as error_message:
            print(f"{error_message}: Files have already been moved")
        
        analyser_path = os.environ["BSS_HOME"] + "analyse_freenrg"
        if not os.path.isfile(f"{unbound}/mbar.txt"):
            try:
                unbound_command = f"{analyser_path} mbar -i {unbound}/lambda*/simfile.dat -o {unbound}/mbar.txt --overlap --subsampling"
                sp.check_output(unbound_command, shell=True)
            
            except sp.CalledProcessError as error_message:
                print(error_message.output)
                print("Trying again without subsampling.")
                warnings.warn(f"Warning: Disabling subsampling may meen results are unreliable. Please check the unbound transformation {transformation}")
                unbound_command = f"{analyser_path} mbar -i {unbound}/lambda*/simfile.dat -o {unbound}/mbar.txt --overlap"

            with open(unbound + "mbar.out", "w") as file:
                sp.run(unbound_command, shell=True, stdout=file)

        if not os.path.isfile(f"{bound}/mbar.txt"):
            try:
                bound_command = f"{analyser_path} mbar -i {bound}/lambda*/simfile.dat -o {bound}/mbar.txt --overlap --subsampling"
                sp.check_output(bound_command, shell=True)
            
            except sp.CalledProcessError as error_message:
                print(error_message.output)
                print("Trying again without subsampling.")
                warnings.warn(f"Warning: Disabling subsampling may meen results are unreliable. Please check the bound transformation {transformation}")
                bound_command = f"{analyser_path} mbar -i {bound}/lambda*/simfile.dat -o {bound}/mbar.txt --overlap"

            with open(bound + "mbar.out", "w") as file:
                sp.run(bound_command, shell=True, stdout=file)

        unbound_free_energy, unbound_error = read_mbar(unbound + "/mbar.txt")
        bound_free_energy, bound_error = read_mbar(bound + "/mbar.txt")
        relative_binding_free_energy = None
        error = None

        try:
            relative_binding_free_energy = bound_free_energy - unbound_free_energy
            error = math.sqrt(bound_error ** 2 + unbound_error ** 2)
        except TypeError as error_message:
            print(f"{transformation}: {error_message}")

        data = [transformation, relative_binding_free_energy, error]
        data_line = ",".join(str(item) for item in data) + "\n"

        data_file = outputs + "/" + engine + f"_{i}_raw.csv"
        with open(data_file, "w+") as output_file:
            lines = output_file.readlines()

 
            if len(lines) == 0:
                header = "transformation,free-energy,error\n"
                output_file.write(header)

        with open(data_file, "r") as input_file:
            input_lines = input_file.readlines()
        
        if data in input_lines:
            warnings.warn(f"Result for {transformation} already in {data_file}")
        
        else:
            with open(data_file, "a") as output_file:
                output_file.write(data_line)


def compute_rmsd(directory, topology_format="PARM7"):
    """
    Compute root mean square deviation in a simulation

    Parameters:
    -----------
    directory: str
        full path to the directory containing lambda directories eg. project/SOMD_1/lig_A_lig_B/bound/
    engine: str
        name of the MD engine used
    topology_format: str
        file format of the topology file, currently only support PARM7 

    Return:
    -------
    time, rmsd_values: tuple
        tuple of two lists containing the time and RMSD values 
    """
    lambda_directories = functions.get_files(directory + "lambda_*/")

    topology_file = functions.get_files(lambda_directories[0] + "/*.prm7")[0]
    
    trajectories = functions.get_files(directory + "/lambda_*/*.dcd") 
    
    with mda.lib.formats.libdcd.DCDFile(trajectories[0]) as trajectory:
        frames = [frame for frame in trajectory]
    first_frame = frames[0].xyz

    reference_universe = mda.Universe(topology_file, first_frame, topology_format=topology_format)
    universe = mda.Universe(topology_file, trajectories, topology_format=topology_format)
    
    reference_coordinates = reference_universe.select_atoms("resname LIG")
    ligand = universe.select_atoms("resname LIG")
    
    rmsd = mda.analysis.rms.RMSD(ligand, reference_coordinates)
    rmsd.run()
    rmsd_result = rmsd.results["rmsd"].T

    time = rmsd_result[0] / 1000
    rmsd_values = rmsd_result[2]

    return time, rmsd_values


def save_rmsds(protocol, transformation):
    
    outputs = functions.path_exists(protocol["outputs"])
    engine = protocol["engine"]
    n_repeats = functions.check_int(protocol["repeats"])

    for repeat in range(1, n_repeats + 1):
        
        path = f"{outputs}/{engine}_{repeat}/{transformation}/"
        unbound = path + "unbound/"
        bound = path + "bound/"    

        unbound_time, unbound_rmsd = compute_rmsd(unbound)
        bound_time, bound_rmsd = compute_rmsd(bound)
        
        unbound_rmsd_array = np.array([unbound_time, unbound_rmsd])
        np.save(path  + "unbound_rmsd.npy", unbound_rmsd_array)
        bound_rmsd_array = np.array([bound_time, bound_rmsd])
        np.save(path  + "bound_rmsd.npy", bound_rmsd_array)


def main():
    
    parser = argparse.ArgumentParser(description="MEZE: MetalloEnZymE FF-builder for alchemistry")

    parser.add_argument("protocol_file",
                        help="protocol file containing equilibration options",
                        type=str,
                        default=os.getcwd() + "/afe/protocol.dat")
    
    parser.add_argument("transformation",
                        help="the pair of ligands undergoing AFE transformation, e.g. ligand_1~ligand_2",
                        type=str)
    
    parser.add_argument("experimental_file",
                        default=os.getcwd() + "/afe/experimental_K_i.csv",
                        help="file containing experimental inhibition constants and errors")
    
    parser.add_argument("-s",
                        "--separator",
                        help="character separating the two ligand names",
                        default="~",
                        type=meze.character)
    
    arguments = parser.parse_args()
    protocol_file = functions.file_exists(arguments.protocol_file)
    protocol = functions.read_protocol(protocol_file)
    transformation = arguments.transformation

    fix_simfile(protocol, transformation)
    save_results(protocol, transformation)
    save_rmsds(protocol, transformation)
    

if __name__ == "__main__":
    main()

