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
import matplotlib.pyplot as plt
import seaborn as sns
logger = logging.getLogger()
logger.setLevel(logging.CRITICAL)


def read_free_energy(mbar_textfile):
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

        unbound_free_energy, unbound_error = read_free_energy(unbound + "/mbar.txt")
        bound_free_energy, bound_error = read_free_energy(bound + "/mbar.txt")
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
        
        input_lines = []
        if os.path.isfile(data_file): 
            with open(data_file, "r") as output_file:
                input_lines = output_file.readlines()

        with open(data_file, "a") as output_file:
            if len(input_lines) == 0:
                header = "transformation,free-energy,error\n"
                output_file.write(header)

        if data_line in input_lines:
            warnings.warn(f"Result for {transformation} already in {data_file}")
        
        else:
            with open(data_file, "a") as output_file:
                print(f"Writing {data_line} to file {output_file}")
                output_file.write(data_line)


def plot_and_compute_rmsd(directory, topology_format="PARM7"):
    """
    Compute root mean square deviation in a simulation

    Parameters:
    -----------
    directory: str
        full path to the directory containing lambda directories e.g. project/SOMD_1/lig_A_lig_B/bound/
    engine: str
        name of the MD engine used
    topology_format: str
        file format of the topology file, currently only support PARM7 

    Return:
    -------
    time, rmsd_values: tuple
        tuple of two lists containing the time and RMSD values 
    """
    topology_files = functions.get_files(directory + "lambda_*/*.prm7") 
    trajectories = functions.get_files(directory + "/lambda_*/*.dcd") 

    for i in range(len(topology_files)):
        save_path, _ = os.path.split(topology_files[i])
        savename = save_path + "/rmsd"
        if not os.path.isfile(savename + ".npy"):
            with mda.lib.formats.libdcd.DCDFile(trajectories[i]) as trajectory:
                frames = [frame for frame in trajectory]
            first_frame = frames[0].xyz

            reference_universe = mda.Universe(topology_files[i], first_frame, topology_format=topology_format)
            universe = mda.Universe(topology_files[i], trajectories[i], topology_format=topology_format)
            
            reference_coordinates = reference_universe.select_atoms("resname LIG")
            ligand = universe.select_atoms("resname LIG")
            
            rmsd = mda.analysis.rms.RMSD(ligand, reference_coordinates)
            rmsd.run()
            rmsd_result = rmsd.results["rmsd"].T

            time = rmsd_result[0] / 1000
            rmsd_values = rmsd_result[2]
                
            save_array = np.array([time, rmsd_values])
            np.save(savename, save_array)           
            
            fig, ax = plt.subplots(1, 1, figsize=(10, 10))
            sns.set(context="notebook", palette="colorblind", style="ticks", font_scale=2)
            ax.plot(time, rmsd_values, "k-")
            ax.set_xlabel("Time (ns)")
            ax.set_ylabel("RMSD ($\AA$)")
            fig.savefig(savename + ".png", dpi=1000)
            plt.close(fig)


def save_rmsds(protocol, transformation):
    """
    Compute RMSD of each trajectory + pairwise lambda RMSDs and save as npy files

    Parameters:
    -----------
    protocol: dict
        protocol file as a dictionary

    transformation: str 
        name of the transformation

    Return:
    -------
    """
    
    outputs = functions.path_exists(protocol["outputs"])
    engine = protocol["engine"]
    n_repeats = functions.check_int(protocol["repeats"])

    for repeat in range(1, n_repeats + 1):
        
        path = f"{outputs}/{engine}_{repeat}/{transformation}/"
        unbound = path + "unbound/"
        bound = path + "bound/"    
        
        plot_and_compute_rmsd(unbound)
        plot_and_compute_rmsd(bound)

        pairwise_unbound_file_first_frame = path + "pairwise_unbound_rmsd_first_frame.npy"
        pairwise_bound_file_first_frame = path + "pairwise_bound_rmsd_first_frame.npy"

        if not os.path.isfile(pairwise_unbound_file_first_frame):
            pairwise_unbound_rmsd_first_frame = compute_pairwise_lambda_rmsd(unbound, frame=0)
            np.save(pairwise_unbound_file_first_frame, pairwise_unbound_rmsd_first_frame)
        
        if not os.path.isfile(pairwise_bound_file_first_frame):
            pairwise_bound_rmsd_first_frame = compute_pairwise_lambda_rmsd(bound, frame=0)
            np.save(pairwise_bound_file_first_frame, pairwise_bound_rmsd_first_frame)

        pairwise_unbound_file_last_frame = path + "pairwise_unbound_rmsd_last_frame.npy"
        pairwise_bound_file_last_frame = path + "pairwise_bound_rmsd_last_frame.npy"

        if not os.path.isfile(pairwise_unbound_file_last_frame):
            pairwise_unbound_rmsd_last_frame = compute_pairwise_lambda_rmsd(unbound, frame=-1)
            np.save(pairwise_unbound_file_last_frame, pairwise_unbound_rmsd_last_frame)
        
        if not os.path.isfile(pairwise_bound_file_last_frame):
            pairwise_bound_rmsd_last_frame = compute_pairwise_lambda_rmsd(bound, frame=-1)
            np.save(pairwise_bound_file_last_frame, pairwise_bound_rmsd_last_frame)



def compute_pairwise_lambda_rmsd(directory, frame):
    """
    Compute the pairwise RMSD between lambda windows at given frame

    Parameters:
    -----------
    directory: str
        full path to the directory containing lambda directories e.g. project/SOMD_1/lig_A_lig_B/bound/

    frame: int
        frame at which pairwise RMSD is calculated

    Return:
    -------
    pairwise_rmsd_matrix: np.array
        symmetric matrix whose dimensions are n_lambdas x n_lambdas; 
        each cell is the pairwise RMSD between two lambda windows calculated at given frame

    """
    frame_index = functions.check_int(frame)

    lambda_directories = functions.get_files(directory + "lambda_*/")

    pairwise_rmsd_matrix = np.zeros(shape=(len(lambda_directories),
                                           len(lambda_directories)))
    
    for i in range(len(lambda_directories)):
        for j in range(len(lambda_directories)):

            lambda_1 = lambda_directories[i]
            lambda_2 = lambda_directories[j]
            topology_1 = functions.get_files(lambda_1 + "*.prm7")[0]
            topology_2 = functions.get_files(lambda_2 + "*.prm7")[0]
            trajectory_1 = functions.get_files(lambda_1 + "*.dcd")[0]
            trajectory_2 = functions.get_files(lambda_2 + "*.dcd")[0]

            universe_1 = mda.Universe(topology_1, trajectory_1, topology_format="PARM7")
            universe_2 = mda.Universe(topology_2, trajectory_2, topology_format="PARM7")

            ligand_1 = universe_1.select_atoms("resname LIG")
            ligand_2 = universe_2.select_atoms("resname LIG")

            universe_1.trajectory[frame_index]
            universe_2.trajectory[frame_index]

            lambda_1 = ligand_1.positions.copy()
            lambda_2 = ligand_2.positions.copy()

            rmsd = MDAnalysis.analysis.rms.rmsd(lambda_1, lambda_2)   

            pairwise_rmsd_matrix[i][j] = rmsd     

    lower_triangle_indices = np.tril_indices(len(pairwise_rmsd_matrix), -1)
    pairwise_rmsd_matrix[lower_triangle_indices] = pairwise_rmsd_matrix.T[lower_triangle_indices]
    
    return pairwise_rmsd_matrix


def save_overlap_matrix(protocol, transformation):
    """
    Read in overlap matrix and save as a numpy array to transformation directory.

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

        unbound_mbar_file = unbound + "/mbar.txt"
        bound_mbar_file = bound + "/mbar.txt"

        unbound_output_file = path + "unbound_overlap_matrix.npy"
        bound_output_file = path + "bound_overlap_matrix.npy"

        if not os.path.isfile(unbound_output_file):
            unbound_overlap_matrix = read_overlap_matrix(unbound_mbar_file)
            np.save(unbound_output_file, unbound_overlap_matrix)
        
        if not os.path.isfile(bound_output_file):
            bound_overlap_matrix = read_overlap_matrix(bound_mbar_file)
            np.save(bound_output_file, bound_overlap_matrix)


def read_overlap_matrix(mbar_file):
    """
    Open mbar.txt and read lines containgin the overlap matrix.
    Return the matrix as a numpy array.

    Parameters:
    -----------
    mbar_file: str
        mbar.txt results file

    Return:
    -------
    matrix: np.array
        overlap matrix as an array
    """
    with open(mbar_file, "r") as file:
        mbar_lines = file.readlines()
    
    start_index = 1
    end_index = -1
    for i in range(len(mbar_lines)):
        if "#Overlap matrix" in mbar_lines[i]:
            start_index = i + 1
        elif "#DG from neighbouring lambda in kcal/mol" in mbar_lines[i]:
            end_index = i
    matrix_lines = mbar_lines[start_index:end_index]        

    matrix_list = []
    for line in matrix_lines:
        split_line = line.replace("\n", "").split()
        new_line = [float(value) for value in split_line]
        matrix_list.append(new_line)
    
    matrix = np.array(matrix_list)
    return matrix


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

    save_results(protocol, transformation)
    save_rmsds(protocol, transformation)    
    save_overlap_matrix(protocol, transformation)


if __name__ == "__main__":
    main()


