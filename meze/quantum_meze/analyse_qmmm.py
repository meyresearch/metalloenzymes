import MDAnalysis as mda
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
import functions
import os
import numpy as np


def get_rmsd(rmsd_file, runtime):
    """
    read in an xvg-formatted file and ouput lists of time, rmsd
    """
    df = pd.read_csv(rmsd_file, sep=r"\s+", header=0)
    frames = df["#Frame"].tolist()
    ns_per_frame = runtime / len(frames)
    time = [frame * ns_per_frame for frame in frames]
    clean_df = df.drop(columns="#Frame")
    clean_df["time"] = time
    return clean_df


def get_bond_lengths_xvg_reader(distance_file):
    """
    get bond lengths between zinc and its coordinating atoms
    """    
    reader = mda.auxiliary.XVG.XVGReader(distance_file)
    time = [step.data[0]/10 for step in reader]
    get_distance = lambda index: [step.data[index] for step in reader]
    distance_1 = get_distance(1)
    distance_2 = get_distance(2)
    distance_3 = get_distance(3)
    distance_4 = get_distance(4)
    return time, distance_1, distance_2, distance_3, distance_4


def get_bond_lengths(distances_file, runtime):

    df = pd.read_csv(distances_file, sep=r"\s+", header=0)
    frames = df["#Frame"].tolist()
    ns_per_frame = runtime / len(frames)
    time = [frame * ns_per_frame for frame in frames]
    clean_df = df.drop(columns="#Frame")
    clean_df["time"] = time
    return clean_df


def get_running_average(time, quantity, window=100):
    """
    calculate the running average over a window of N ps
    """
    dataframe = pd.DataFrame({"time": time, "distance": quantity})
    return dataframe["distance"].rolling(window=window).mean().to_numpy()


def main():

    parser = argparse.ArgumentParser(description="solvation for meze workflow")

    parser.add_argument("-bb",
                        "--backbone-rmsd-file",
                        dest="backbone_file",
                        type=str,
                        help="datafile containing protein backbone rmsd")
    
    parser.add_argument("-m",
                        "--metal-site-rmsd-file",
                        dest="metal_site_file",
                        type=str,
                        help="datafile containing metal site residues' rmsd")

    arguments = parser.parse_args()

    metal_rmsd_file = functions.file_exists(arguments.metal_site_file)

    rmsd_df = get_rmsd(metal_rmsd_file, runtime=10)


    fig, ax = plt.subplots(figsize=(10, 10))
    sns.set_theme(context="notebook", palette="colorblind", style="ticks", font="DejaVu Sans")

    ax.plot(rmsd_df["time"], rmsd_df["h85"], c="#D0006F", label="H85")
    ax.plot(rmsd_df["time"], rmsd_df["h87"], c="#0099AB", label="H87")
    ax.plot(rmsd_df["time"], rmsd_df["h150"], c="#29C2DE", label="H150")
    ax.plot(rmsd_df["time"], rmsd_df["d89"], c="#830065", label="D89")
    ax.plot(rmsd_df["time"], rmsd_df["h169"], c="#D0006F", label="C169")
    ax.plot(rmsd_df["time"], rmsd_df["h211"], c="#C25E03", label="H211")
    ax.plot(rmsd_df["time"], rmsd_df["h239"], c="#F9A800", label="Lig 1")
    ax.set_xlabel("time (ps)", fontsize=14)
    ax.set_ylabel("RMSD (Å)", fontsize=14)
    plt.legend(frameon=False, fontsize=14)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    # sns.despine()
    plt.savefig("rmsd_metal_site.png", dpi=800)


if __name__ == "__main__":
    main()