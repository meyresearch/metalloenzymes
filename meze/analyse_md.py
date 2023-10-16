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
    reader = mda.auxiliary.XVG.XVGReader(rmsd_file)

    frames = [step.data[0]  for step in reader]
    ns_per_frame = runtime / len(frames)
    time = [frame * ns_per_frame for frame in frames]
    rmsd = [step.data[1] for step in reader]
    return time, rmsd


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
    
    parser.add_argument("-d",
                        "--distances-file",
                        dest="distances_file",
                        type=str,
                        help="datafile containing metal - residue distances")
    
    parser.add_argument("protocol_file",
                        help="protocol file containing equilibration options",
                        type=str,
                        default=os.getcwd() + "/qmmm_input_files/protocol.dat")

    arguments = parser.parse_args()
    protocol = functions.input_to_dict(arguments.protocol_file)

    bb_rmsd_file = functions.file_exists(arguments.backbone_file)
    metal_rmsd_file = functions.file_exists(arguments.metal_site_file)

    bb_time, bb_rmsd = get_rmsd(bb_rmsd_file, protocol["sampling time"])
    metal_time, metal_rmsd = get_rmsd(metal_rmsd_file, protocol["sampling time"])

    fig, ax = plt.subplots(figsize=(10, 10))
    sns.set_theme(context="notebook", palette="colorblind", style="ticks", font="DejaVu Sans")
    ax.plot(bb_time, bb_rmsd, c="#44AA99", label="Backbone")
    ax.plot(metal_time, metal_rmsd, c="#999933", label="Metal ligands")
    ax.set_xlabel("time (ns)", fontsize=14)
    ax.set_ylabel("RMSD (Å)", fontsize=14)
    plt.legend(frameon=False, fontsize=14)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    # sns.despine()
    plt.savefig("rmsd_backbone_metal_site.png", dpi=800)

    distances_file = functions.file_exists(arguments.distances_file)
    distances_df = get_bond_lengths(distances_file, protocol["sampling time"])
    # distances = distances_df.to_dict(orient="list")

    # avg = np.mean(distances_df["lig_o1"])
    # print(avg)
    # std = np.std(distances_df["lig_o2"])
    # print(std)

    means = distances_df.mean()
    stds = distances_df.std()
    print(f"MEANS: \n{means}")
    print(f"STDS: \n{stds}")

    sns.set(context="notebook", palette="colorblind", style="ticks", font_scale=2)

    plot = sns.jointplot(data=distances_df, x="time", y="lig_o1", kind="scatter", color="#0099AB",)
    plot.ax_joint.cla()
    plt.sca(plot.ax_joint)
    plt.scatter(distances_df["time"], distances_df["lig_o1"], s=0)
    plot.ax_joint.plot(distances_df["time"], distances_df["lig_o1"], color="black", alpha=0.8)

    plt.gcf().set_size_inches(10, 10)
    plot.set_axis_labels("Time (ns)", "Distance (Å)")
    plt.setp(plot.ax_marg_x.patches, color="w")
    # plt.setp(plot.ax_marg_x.spines, color="w")
    plot.ax_marg_x.remove() 
    plot.ax_joint.spines["top"].set_visible(True)
    sns.despine()
    plt.tight_layout()
    plt.savefig("lig_o1.png", dpi=800)



if __name__ == "__main__":
    main()


