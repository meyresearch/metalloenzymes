import matplotlib.pyplot as plt
import MDAnalysis as mda
import pandas as pd
import numpy as np
import argparse
import seaborn as sns


def get_rmsd(rmsd_file: str) -> tuple:
    """
    read in an xvg-formatted file and ouput lists of time, rmsd
    """
    reader = mda.auxiliary.XVG.XVGReader(rmsd_file)
    time = [step.data[0] for step in reader]
    rmsd = [step.data[1] for step in reader]
    return time, rmsd


def get_bond_lengths(distance_file: str) -> tuple:
    """
    get bond lengths between zinc and its coordinating atoms
    """    
    reader = mda.auxiliary.XVG.XVGReader(distance_file)
    time = [step.data[0] for step in reader]
    get_distance = lambda index: [step.data[index] for step in reader]
    distance_1 = get_distance(1)
    distance_2 = get_distance(2)
    distance_3 = get_distance(3)
    distance_4 = get_distance(4)
    return time, distance_1, distance_2, distance_3, distance_4


def get_running_average(time: list, quantity: list, window=100) -> np.array:
    """
    calculate the running average over a window of N ps
    """
    dataframe = pd.DataFrame({"time": time, "distance": quantity})
    return dataframe["distance"].rolling(window=window).mean().to_numpy()


def plot_distance(time, bond_distance, ligand_number, distance_name, savedir):

    average_distance = np.mean(bond_distance)
    std_distance = np.std(bond_distance)
    running_average = get_running_average(time, bond_distance, window=1000)
    distance_file = f"{savedir}/ligand_{ligand_number}_{distance_name}.txt"
    
    with open(distance_file, "a") as file:
        file.write(rf"{distance_name}: {average_distance:.4f} +/- {std_distance:.4} $\AA$")

    data_array = np.array([time, bond_distance])
    data_array_tr = np.transpose(data_array)
    dataframe = pd.DataFrame(data_array_tr, columns=["time", "coordination"])
    
    sns.set(context="notebook", palette="colorblind", style="ticks", font_scale=2)
    plot = sns.jointplot(data=dataframe, x="time", y="coordination", kind="scatter", color="#0099AB")
    plot.ax_joint.cla()
    plt.sca(plot.ax_joint)
    plt.scatter(dataframe["time"], dataframe["coordination"], s=0)
    plot.ax_joint.plot(dataframe["time"], dataframe["coordination"], color="gray", alpha=0.7)
    plt.setp(plot.ax_marg_x.patches, color="w")

    # if distance_name == "zn2_d":
    #     plt.hlines(1.99, time[0], time[-1], linestyles="--", colors="black", label="Reference")
    #     # plt.fill_between(time, 1.99 + 0.5, 1.99 - 0.5, color="grey", alpha=0.5)
    # if distance_name == "zn2_c":
    #     plt.hlines(2.31, time[0], time[-1], linestyles="--", colors="black", label="Reference")
    #     # plt.fill_between(time, 2.31 + 0.5, 2.31 - 0.5, color="grey", alpha=0.5)
    # if distance_name == "zn2_h":
    #     plt.hlines(2.03, time[0], time[-1], linestyles="--", colors="black", label="Reference")
    #     # plt.fill_between(time, 2.03 + 0.5, 2.03 - 0.5, color="grey", alpha=0.5)
    # if distance_name == "zn2_lig":
    #     plt.hlines(2.09, time[0], time[-1], linestyles="--", colors="black", label="Reference")
    #     # plt.fill_between(time, 2.09 + 0.5, 2.09 - 0.5, color="grey", alpha=0.5)
    
    plt.plot(time, running_average, color="#0099AB")
    plot.ax_marg_x.remove() 
    plot.ax_joint.spines["top"].set_visible(True)
    plt.gcf().set_size_inches(8, 8)
    plot.set_axis_labels("Time (ns)", "Distance (Å)")
    sns.despine()
    # plt.legend(loc="upper right")
    plt.tight_layout()
    plt.savefig(f"{savedir}/ligand_{ligand_number}_{distance_name}_tr.png", transparent=True)

    return plt.gca()


parser = argparse.ArgumentParser(description="analyse ligand MD")

parser.add_argument("-n", "--ligand-number",
                    dest="ligand_number", 
                    type=str,
                    help="ligand number for analysis")

arguments = parser.parse_args()

try:
    n = int(arguments.ligand_number)
except ValueError as e:
    print(f"{e}: This is not a valid ligand number")

distance_file = f"/home/jguven/projects/metalloenzymes/nonbonded_model_vim2/ligand_{n}/production/dist_md.txt"

data = pd.read_csv(distance_file, header=None, sep="\s+", names=["t", "h1", "h2", "h3", "w", "d", "c", "h", "o"])
time = data["t"].to_numpy() * (0.002/1000)

get_distance = lambda column: data[column].to_numpy()

zn1_h1, zn1_h2, zn1_h3, zn1_wat = get_distance("h1"), get_distance("h2"), get_distance("h3"), get_distance("w")

zn2_d, zn2_c, zn2_h, zn2_lig = get_distance("d"), get_distance("c"), get_distance("h"), get_distance("o")

plots_dir = "/home/jguven/projects/metalloenzymes/nonbonded_model_vim2/plots/"

# print("HISTIDINE HID 84 ")

plot1 = plot_distance(time, zn1_h1, n, "zn1_h1", plots_dir)
plot2 = plot_distance(time, zn1_h2, n, "zn1_h2", plots_dir)
plot3 = plot_distance(time, zn1_h3, n, "zn1_h3", plots_dir)
plot4 = plot_distance(time, zn1_wat, n, "zn1_wat", plots_dir)
plot5 = plot_distance(time, zn2_d, n, "zn2_d", plots_dir)
plot6 = plot_distance(time, zn2_c, n, "zn2_c", plots_dir)
plot7 = plot_distance(time, zn2_h, n, "zn2_h", plots_dir)
plot8 = plot_distance(time, zn2_lig, n, "zn2_lig", plots_dir)


# h1_avg = np.mean(zn1_h1)
# h1_std = np.std(zn1_h1)
# runnning_avg_h1 = get_running_average(time, zn1_h1, window=1000)
# print(f"{h1_avg} +/- {h1_std}")
# data1 = np.array([time, zn1_h1])
# data1_tr = np.transpose(data1)
# df1 = pd.DataFrame(data1_tr, columns=["time", "ligand"])
# sns.set(context="notebook", palette="colorblind", style="ticks", font_scale=2)
# plot = sns.jointplot(data=df1, x="time", y="ligand", kind="scatter", color="#0099AB",)
# plot.ax_joint.cla()
# plt.sca(plot.ax_joint)
# plt.scatter(df1["time"], df1["ligand"], s=0)
# plot.ax_joint.plot(df1["time"], df1["ligand"], color="gray", alpha=0.5)
# plt.plot(time, runnning_avg_h1, color="#0099AB")
# plt.setp(plot.ax_marg_x.patches, color="w")
# plot.ax_marg_x.remove() 
# plot.ax_joint.spines["top"].set_visible(True)
# plt.gcf().set_size_inches(8, 8)
# plot.set_axis_labels("Time (ns)", "Distance (Å)")
# sns.despine()
# plt.tight_layout()
# plt.savefig(f"{plots_dir}/ligand_{n}_zn1_hid84.pdf")

# print("HISTIDINE HIE 86 ")
# h2_avg = np.mean(zn1_h2)
# h2_std = np.std(zn1_h2)
# runnning_avg_h2 = get_running_average(time, zn1_h2, window=1000)
# print(f"{h2_avg} +/- {h2_std}")
# data2 = np.array([time, zn1_h2])
# data2_tr = np.transpose(data2)
# df2 = pd.DataFrame(data2_tr, columns=["time", "ligand"])
# sns.set(context="notebook", palette="colorblind", style="ticks", font_scale=2)
# plot = sns.jointplot(data=df2, x="time", y="ligand", kind="scatter", color="#0099AB",)
# plot.ax_joint.cla()
# plt.sca(plot.ax_joint)
# plt.scatter(df2["time"], df2["ligand"], s=0)
# plot.ax_joint.plot(df2["time"], df2["ligand"], color="gray", alpha=0.5)
# plt.plot(time, runnning_avg_h2, color="#0099AB")
# plt.setp(plot.ax_marg_x.patches, color="w")
# plot.ax_marg_x.remove() 
# plot.ax_joint.spines["top"].set_visible(True)
# plt.gcf().set_size_inches(8, 8)
# plot.set_axis_labels("Time (ns)", "Distance (Å)")
# sns.despine()
# plt.tight_layout()
# plt.savefig(f"{plots_dir}/ligand_{n}_zn1_hie86.pdf")

# print("HISTIDINE HID 149 ")
# h3_avg = np.mean(zn1_h3)
# h3_std = np.std(zn1_h3)
# runnning_avg_h3 = get_running_average(time, zn1_h3, window=1000)
# print(f"{h3_avg} +/- {h3_std}")
# data3 = np.array([time, zn1_h3])
# data3_tr = np.transpose(data3)
# df3 = pd.DataFrame(data3_tr, columns=["time", "ligand"])
# sns.set(context="notebook", palette="colorblind", style="ticks", font_scale=2)
# plot = sns.jointplot(data=df3, x="time", y="ligand", kind="scatter", color="#0099AB",)
# plot.ax_joint.cla()                                                                                                                             
# plt.sca(plot.ax_joint)
# plt.scatter(df3["time"], df3["ligand"], s=0)
# plot.ax_joint.plot(df3["time"], df3["ligand"], color="gray", alpha=0.5)
# plt.plot(time, runnning_avg_h3, color="#0099AB")
# plt.setp(plot.ax_marg_x.patches, color="w")
# plot.ax_marg_x.remove() 
# plot.ax_joint.spines["top"].set_visible(True)
# plt.gcf().set_size_inches(8, 8)
# plot.set_axis_labels("Time (ns)", "Distance (Å)")
# sns.despine()
# plt.tight_layout()
# plt.savefig(f"{plots_dir}/ligand_{n}_zn1_hid149.pdf")

# print("XTAL WATER")
# wat_avg = np.mean(zn1_wat)
# wat_std = np.std(zn1_wat)
# running_avg_wat = get_running_average(time, zn1_wat, window=1000)
# print(f"{wat_avg} +/- {wat_std}")
# data4 = np.array([time, zn1_wat])
# data4_tr = np.transpose(data4)
# df4 = pd.DataFrame(data4_tr, columns=["time", "ligand"])
# sns.set(context="notebook", palette="colorblind", style="ticks", font_scale=2)
# plot = sns.jointplot(data=df4, x="time", y="ligand", kind="scatter", color="#0099AB",)
# plot.ax_joint.cla()
# plt.sca(plot.ax_joint)
# plt.scatter(df4["time"], df4["ligand"], s=0)
# plot.ax_joint.plot(df4["time"], df4["ligand"], color="gray", alpha=0.5)
# plt.plot(time, running_avg_wat, color="#0099AB")
# plt.setp(plot.ax_marg_x.patches, color="w")
# plot.ax_marg_x.remove() 
# plot.ax_joint.spines["top"].set_visible(True)
# plt.gcf().set_size_inches(8, 8)
# plot.set_axis_labels("Time (ns)", "Distance (Å)")
# sns.despine()
# plt.tight_layout()
# plt.savefig(f"{plots_dir}/ligand_{n}_zn1_wat.pdf")

# print("ASPARTIC ACID ASP 88 ")
# d_avg = np.mean(zn2_d)
# d_std = np.std(zn2_d)
# running_avg_d = get_running_average(time, zn2_d, window=1000)
# print(f"{d_avg} +/- {d_std}")
# data5 = np.array([time, zn2_d])
# data5_tr = np.transpose(data5)
# df5 = pd.DataFrame(data5_tr, columns=["time", "ligand"])
# sns.set(context="notebook", palette="colorblind", style="ticks", font_scale=2)
# plot = sns.jointplot(data=df5, x="time", y="ligand", kind="scatter", color="#0099AB",)
# plot.ax_joint.cla()
# plt.sca(plot.ax_joint)
# plt.scatter(df5["time"], df5["ligand"], s=0)
# plot.ax_joint.plot(df5["time"], df5["ligand"], color="gray", alpha=0.5)
# plt.plot(time, running_avg_d, color="#0099AB")
# plt.setp(plot.ax_marg_x.patches, color="w")
# plot.ax_marg_x.remove() 
# plot.ax_joint.spines["top"].set_visible(True)
# plt.gcf().set_size_inches(8, 8)
# plot.set_axis_labels("Time (ns)", "Distance (Å)")
# sns.despine()
# plt.tight_layout()
# plt.savefig(f"{plots_dir}/ligand_{n}_zn2_asp88.pdf")

# print("CYSTEINE 168")
# c_avg = np.mean(zn2_c)
# c_std = np.std(zn2_c)
# running_avg_c = get_running_average(time, zn2_c, window=1000)
# print(f"{c_avg} +/- {c_std}")
# data6 = np.array([time, zn2_c])
# data6_tr = np.transpose(data6)
# df6 = pd.DataFrame(data6_tr, columns=["time", "ligand"])
# sns.set(context="notebook", palette="colorblind", style="ticks", font_scale=2)
# plot = sns.jointplot(data=df6, x="time", y="ligand", kind="scatter", color="#0099AB",)
# plot.ax_joint.cla()
# plt.sca(plot.ax_joint)
# plt.scatter(df6["time"], df6["ligand"], s=0)
# plot.ax_joint.plot(df6["time"], df6["ligand"], color="gray", alpha=0.5)
# plt.plot(time, running_avg_c, color="#0099AB")
# plt.setp(plot.ax_marg_x.patches, color="w")
# plot.ax_marg_x.remove() 
# plot.ax_joint.spines["top"].set_visible(True)
# plt.gcf().set_size_inches(8, 8)
# plot.set_axis_labels("Time (ns)", "Distance (Å)")
# sns.despine()
# plt.tight_layout()
# plt.savefig(f"{plots_dir}/ligand_{n}_zn2_cys168.pdf")

# print("HISTIDINE 210")
# h_avg = np.mean(zn2_h)
# h_std = np.std(zn2_h)
# running_avg_h = get_running_average(time, zn2_h, window=1000)
# print(f"{h_avg} +/- {h_std}")
# data7 = np.array([time, zn2_h])
# data7_tr = np.transpose(data7)
# df7 = pd.DataFrame(data7_tr, columns=["time", "ligand"])
# sns.set(context="notebook", palette="colorblind", style="ticks", font_scale=2)
# plot = sns.jointplot(data=df7, x="time", y="ligand", kind="scatter", color="#0099AB",)
# plot.ax_joint.cla()
# plt.sca(plot.ax_joint)
# plt.scatter(df7["time"], df7["ligand"], s=0)
# plot.ax_joint.plot(df7["time"], df7["ligand"], color="gray", alpha=0.5)
# plt.plot(time, running_avg_h, color="#0099AB")
# plt.setp(plot.ax_marg_x.patches, color="w")
# plot.ax_marg_x.remove() 
# plot.ax_joint.spines["top"].set_visible(True)
# plt.gcf().set_size_inches(8, 8)
# plot.set_axis_labels("Time (ns)", "Distance (Å)")
# sns.despine()
# plt.tight_layout()
# plt.savefig(f"{plots_dir}/ligand_{n}_zn2_hid210.pdf")

# print("LIGAND")
# lig_avg = np.mean(zn2_lig)
# lig_std = np.std(zn2_lig)
# running_avg_lig = 
# print(f"{lig_avg} +/- {lig_std}")
# data8 = np.array([time, zn2_lig])
# data8_tr = np.transpose(data8)
# df8 = pd.DataFrame(data8_tr, columns=["time", "ligand"])
# sns.set(context="notebook", palette="colorblind", style="ticks", font_scale=2)
# plot = sns.jointplot(data=df8, x="time", y="ligand", kind="scatter", color="#0099AB",)
# plot.ax_joint.cla()
# plt.sca(plot.ax_joint)
# plt.scatter(df8["time"], df8["ligand"], s=0)
# plot.ax_joint.plot(df8["time"], df8["ligand"], color="black", alpha=0.8)
# plt.setp(plot.ax_marg_x.patches, color="w")
# plot.ax_marg_x.remove() 
# plot.ax_joint.spines["top"].set_visible(True)
# plt.gcf().set_size_inches(8, 8)
# plot.set_axis_labels("Time (ns)", "Distance (Å)")
# sns.despine()
# plt.tight_layout()
# plt.savefig(f"{plots_dir}/ligand_{n}_zn2_lig.pdf")

