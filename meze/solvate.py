"""
Solvate bound and unbound legs of AFE workflow
"""
import BioSimSpace as bss


def solvate_meze(idx, Protein, Network, AFE):
    """
    Solvate unbound and bound systems.

    Parameters:
    -----------
    idx: int
        Ligand index for sorting through Network.names
    Protein: 
        Protein class object
    Network: 
        Network class object
    AFE: 
        AlchemicalFreeEnergy class object

    Return:
    -------
    """
    ligands = Network.ligands
    names = Network.names
    ligand_number = Network.names[idx].split("_")[-1]
    print(f"Solvating unbound ligand {ligand_number}")
    ligand_parameters = ligands[idx].parameterise(Network.forcefield, Network.charge)
    unbound_box, unbound_box_angles = AFE.create_box(ligand_parameters)
    solvated_ligand = bss.Solvent.solvate(model=Protein.water_model, 
                                            molecule=ligand_parameters, 
                                            box=unbound_box,
                                            angles=unbound_box_angles)
    print(f"Solvating bound ligand {ligand_number}")        
    system_parameters = ligand_parameters + Protein.get_prepared_protein()
    bound_box, bound_box_angles = AFE.create_box(system_parameters)
    solvated_system = bss.Solvent.solvate(model=Protein.water_model,
                                            molecule=system_parameters,
                                            box=bound_box,
                                            angles=bound_box_angles)
    ligand_savename = Network.path + "ligand_" + ligand_number + "_solvated"
    system_savename = Protein.path + "system_" + ligand_number + "_solvated"
    bss.IO.saveMolecules(ligand_savename, solvated_ligand, ["PRM7", "RST7"])
    bss.IO.saveMolecules(system_savename, solvated_system, ["PRM7", "RST7"])
