import cinnabar
import cinnabar.plotting
import pandas as pd
import scipy 
import numpy as np

home_dir = "/backups/jguven/kpc2_most_recent_results/old_config_partially_protonated/"
afe_dir = f"{home_dir}/afe/"
outputs_dir = f"{home_dir}/outputs/"

exp_file = f"{afe_dir}/experimental_K_i.csv"
calculated_file = f"{outputs_dir}/SOMD_results.csv"


BOLTZMANN_CONSTANT = scipy.constants.Boltzmann
AVOGADROS_NUMBER = scipy.constants.Avogadro

def ki_to_dg(ki):
    return (BOLTZMANN_CONSTANT * 300 * AVOGADROS_NUMBER / 4184) * np.log(ki / 1) 


def propagate_exp_ki_error(ki_error, ki):
    return (BOLTZMANN_CONSTANT * 300 * AVOGADROS_NUMBER / 4184) * (ki_error / ki)

exp_df = pd.read_csv(exp_file)

ki_array = exp_df["K_i"].to_numpy()
ki_err_array = exp_df["K_i_err"].to_numpy()

exp_dg = np.array([ki_to_dg(ki) for ki in ki_array])
exp_dg_error = np.array([propagate_exp_ki_error(ki_err_array[i], ki_array[i]) for i in range(len(ki_array))])

ligand_numbers = exp_df["ligand"].tolist()
ligands = [f"ligand_{n}" for n in ligand_numbers]

calculated_df = pd.read_csv(calculated_file, skiprows=10, )

transformations = calculated_df["transformation"].tolist()
calculated_ddg = calculated_df["average_ddg"].to_numpy()
calculated_error = calculated_df["standard_deviation"].to_numpy()
split_transformations = [trans.replace("~", ",") for trans in transformations]

exp_header = ["# Experimental block\n",
              "# Ligand, expt_DG, expt_dDG\n"]
calculated_header = ["# Calculated block\n",
                     "# Ligand1,Ligand2, calc_DDG, calc_dDDG, calc_dDDG(additional)\n"]

# with open("cinnabar_test_kpc2.csv", "w") as file:
#     file.writelines(exp_header)
#     for i in range(len(exp_dg)):
#         exp_data = [ligands[i], exp_dg[i], exp_dg_error[i]]
#         exp_data_line = ", ".join(str(item) for item in exp_data) + "\n"
#         file.write(exp_data_line)



# with open("cinnabar_test_kpc2.csv", "a") as file:
#     file.write("\n")
#     file.writelines(calculated_header)

#     for i in range(len(calculated_ddg)):
#         calculated_data = [split_transformations[i], calculated_ddg[i], calculated_error[i], 0.0]
#         calculated_data_line = ", ".join(str(item) for item in calculated_data) + "\n"
#         file.write(calculated_data_line)


free_energy_map = cinnabar.FEMap.from_csv("/home/jguven/projects/alchemistry/vim2/deprotonated_ligand/outputs/cinnabar_test_vim2.csv")


fe_graph = free_energy_map.to_legacy_graph()

for node, edge in zip(fe_graph.nodes(), fe_graph.edges()):
    print(f"{node}: {edge}")

cinnabar.plotting.plot_DGs(fe_graph, target_name="VIM-2", title="Test")

