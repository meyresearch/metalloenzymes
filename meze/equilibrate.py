"""
Minimise and equilibrate bound and unbound stages.
"""
from definitions import PICOSECOND, KELVIN, ATM
import functions
import BioSimSpace as bss
import AlchemicalFreeEnergy as afe



def unbound(idx, Network, AFE):


    ligand_number = Network.names[idx].split("_")[-1]
    solvated_ligand = Network.ligands[idx]
    
    unbound_directory = functions.mkdir(AFE.path + f"/equilibration/unbound/ligand_{ligand_number}/")
    

    equilibration = afe.EquilibrationSimulation(path=unbound_directory)
    minimised_system = equilibration.minimise(system=solvated_ligand)

    equilibration.set_name("r_nvt")
    equilibration.set_restraints("all")
    equilibration.set_configuration(["dt = 0.0005"])
    restrained_nvt = equilibration.equilibrate(system=minimised_system)
    




def heat_meze_serial(idx, Protein, Network, AFE):
    """
    Run minimisation and equilibration in unbound and bound stages

    Parameters:
    -----------
    idx: int
        ligand index
    Protein: Protein
        Protein class object
    Network: Network
        Network class object
    AFE: AlchemicalFreeEnergy
        AFE Class object

    Return:
    -------
    """

    ligand_number = Network.names[idx].split("_")[-1]
    unbound_directory = AFE.create_directory(AFE.equilibration_dir + f"/unbound/ligand_{ligand_number}/")
    solvated_ligand = Network.ligands[idx]

    unbound_directories = lambda step: functions.mkdir(unbound_directory+step)
    unbound_min = unbound_directories("min")
    unbound_r_nvt = unbound_directories("r_nvt")
    unbound_nvt = unbound_directories("nvt")
    unbound_r_npt = unbound_directories("r_npt")
    unbound_npt = unbound_directories("npt")    

    minimised_ligand = minimise(system=solvated_ligand, 
                                workdir=unbound_min, 
                                AFE=AFE)
    unbound_restrained_nvt = equilibrate(system=minimised_ligand, 
                                         workdir=unbound_r_nvt, 
                                         name="r_nvt",
                                         runtime=AFE.short_nvt,
                                         config=["dt = 0.0005"],
                                         start_t=0, end_t=300, 
                                         restraints="all")
    unbound_nvt = equilibrate(system=unbound_restrained_nvt,
                              workdir=unbound_nvt,
                              name="nvt",
                              runtime=AFE.nvt,
                              temperature=300)
    unbound_restrained_npt = equilibrate(system=unbound_nvt,
                                         workdir=unbound_r_npt,
                                         name="r_npt",
                                         runtime=AFE.npt,
                                         pressure=1, temperature=300,
                                         restraints="heavy")
    equilibrated_ligand = equilibrate(system=unbound_restrained_npt,
                                      workdir=unbound_npt,
                                      name="npt",
                                      runtime=AFE.npt,
                                      pressure=1, temperature=300)
    unbound_savename = unbound_npt + f"/ligand_{ligand_number}"
    bss.IO.saveMolecules(filebase=unbound_savename, system=equilibrated_ligand, fileformat=["PRM7", "RST7"]) 
    Network.ligands[idx] = equilibrated_ligand

    bound_directory = AFE.create_directory(AFE.equilibration_dir + f"/bound/ligand_{ligand_number}/")
    solvated_system = Protein.complex
    bound_directories = lambda step: functions.mkdir(bound_directory+step)
    bound_min = bound_directories("min")
    bound_r_nvt = bound_directories("r_nvt")
    bound_bb_r_nvt = bound_directories("bb_r_nvt")
    bound_nvt = bound_directories("nvt")
    bound_r_npt = bound_directories("r_npt")
    bound_npt = bound_directories("npt") 

    minimised_system = minimise(system=solvated_system, workdir=bound_min, AFE=AFE)
    bound_restrained_nvt = equilibrate(system=minimised_system,
                                       workdir=bound_r_nvt,
                                       name="r_nvt",
                                       runtime=AFE.short_nvt,
                                       config=["dt = 0.0005"],
                                       start_t=0, end_t=300,
                                       restraints="all")
    bound_backbone_restrained_nvt = equilibrate(system=bound_restrained_nvt,
                                                workdir=bound_bb_r_nvt,
                                                name="bb_r_nvt",
                                                runtime=AFE.nvt,
                                                temperature=300,
                                                restraints="backbone")
    bound_nvt_equilibration = equilibrate(system=bound_backbone_restrained_nvt,
                                          workdir=bound_nvt,
                                          name="nvt",
                                          runtime=AFE.nvt,
                                          temperature=300)
    bound_restrained_npt = equilibrate(system=bound_nvt_equilibration,
                                       workdir=bound_r_npt,
                                       name="r_npt",
                                       runtime=AFE.npt,
                                       pressure=1, temperature=300,
                                       restraints="heavy")
    equilibrated_system = equilibrate(system=bound_restrained_npt,
                                      workdir=bound_npt,
                                      name="npt",
                                      runtime=AFE.npt,
                                      pressure=1, temperature=300)
    bound_savename = bound_npt + f"/system_{ligand_number}"
    bss.IO.saveMolecules(filebase=bound_savename, system=equilibrated_system, fileformat=["PRM7", "RST7"])
    Network.bound[idx] = equilibrated_system

    return Network

def main():
    pass


if __name__ == "__main__":
    main()