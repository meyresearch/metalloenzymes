import argparse
import BioSimSpace as bss
import Ligand
from definitions import ANGSTROM
import Network
import functions
import os


def solvate_unbound(network, ligand_name):
    """
    Solvate unbound systems.

    Parameters:
    -----------
    network: Network
        a prepared Network object   
    ligand_name: str
        ligand name
    Return:
    -------
    solvated_ligand: Ligand
        (solvated) Ligand object whose file attribute is the prm7 and rst7 files 
    """ 
    ligand = network.get_ligand_by_name(ligand_name)
    print(f"Solvating unbound ligand {ligand_name}")
    ligand_parameters = ligand.parameterise(network.ligand_forcefield, network.ligand_charge)
    unbound_box, unbound_box_angles = network.create_box(ligand_parameters)
    solvated_molecule = bss.Solvent.solvate(model=network.protein.water_model, 
                                            molecule=ligand_parameters, 
                                            box=unbound_box,
                                            angles=unbound_box_angles)
    ligand_savename = network.ligand_path + ligand_name + "_solvated"
    solvated_files = bss.IO.saveMolecules(ligand_savename, solvated_molecule, ["PRM7", "RST7"])
    solvated_ligand = Ligand.Ligand(file=solvated_files, parameterised=True)
    return solvated_ligand


def solvate_bound(network, ligand_name):
    """
    Solvate bound systems.

    Parameters:
    -----------
    network: Network
        a prepared Network object    
    ligand_name: str
        name of the ligand
    Return:
    -----
    solvated_system: Ligand
        (solvated) Ligand object whose file attribute is the prm7 and rst7 files         
    """
    ligand = network.get_ligand_by_name(ligand_name)
    print(f"Solvating bound ligand {ligand_name}")
    ligand_parameters = ligand.parameterise(network.ligand_forcefield, network.ligand_charge)    
    protein = network.protein.get_molecule()
    system_parameters = ligand_parameters + protein
    bound_box, bound_box_angles = network.create_box(system_parameters)
    solvated_molecules = bss.Solvent.solvate(model=network.protein.water_model,
                                                molecule=system_parameters,
                                                box=bound_box,
                                                angles=bound_box_angles)
    
    system_savename = network.protein_path + "bound_" + ligand_name + "_solvated"
    solvated_files = bss.IO.saveMolecules(system_savename, solvated_molecules, ["PRM7", "RST7"])
    solvated_system = Ligand.Ligand(file=solvated_files, parameterised=True)
    return solvated_system


def create_box(network, molecule):
    """
    Create a bss.Box object for solvation.

    Parameters:
    -----------
    network: Network
        a prepared Network object
    molecule: bss.Molecule
        usually either a protein or a ligand

    Return:
    -------
    tuple: 
        bss.Box and angles
    """

    box_min, box_max = molecule.getAxisAlignedBoundingBox()
    box_size = [y - x for x, y in zip(box_min, box_max)]
    box_area = [x + int(network.box_edges) * ANGSTROM for x in box_size]
    network.box, network.box_angles = None, None
    if network.box_shape == "cubic":
        network.box, network.box_angles = bss.Box.cubic(max(box_area))
    elif network.box_shape == "rhombicDodecahedronHexagon":
        network.box, network.box_angles = bss.Box.rhombicDodecahedronHexagon(max(box_area))
    elif network.box_shape == "rhombicDodecahedronSquare":
        network.box, network.box_angles = bss.Box.rhombicDodecahedronSquare(max(box_area))
    elif network.box_shape == "truncatedOctahedron":
        network.box, network.box_angles = bss.Box.truncatedOctahedron(max(box_area))
    else:
        print(f"Box shape {network.box_shape} not supported.")
    return network.box, network.box_angles


def main():

    parser = argparse.ArgumentParser(description="solvation for meze workflow")

    parser.add_argument("ligand_name",
                        help="ligand name",
                        type=str)
    
    parser.add_argument("protocol_file",
                        help="protocol file containing equilibration options",
                        type=str,
                        default=os.getcwd() + "/afe/protocol.dat")
    
    arguments = parser.parse_args()
    
    protocol = functions.input_to_dict(arguments.protocol_file)

    network = Network.Network(prepared=True,
                              workdir=protocol["project directory"],
                              ligand_path=protocol["ligand directory"],
                              group_name=protocol["group name"],
                              protein_file=protocol["prepared protein file"],
                              protein_path=protocol["protein directory"],
                              outputs=protocol["outputs"],
                              afe_input_path=protocol["afe input directory"],
                              equilibration_path=protocol["equilibration"],
                              water_model=protocol["water model"],
                              ligand_ff=protocol["ligand forcefield"],
                              protein_ff=protocol["protein forcefield"],
                              ligand_charge=protocol["ligand charge"],
                              engine=protocol["engine"],
                              sampling_time=protocol["sampling time"],
                              box_edges=protocol["box edges"],
                              box_shape=protocol["box shape"],
                              min_steps=protocol["minimisation steps"],
                              short_nvt=protocol["short nvt"],
                              nvt=protocol["nvt"],
                              npt=protocol["npt"],
                              min_dt=protocol["minimisation stepsize"],
                              min_tol=protocol["minimisation tolerance"],
                              repeats=protocol["repeats"],
                              temperature=protocol["temperature"],
                              pressure=protocol["pressure"])

    solvate_unbound(network=network, ligand_name=arguments.ligand_name)
    solvate_bound(network=network, ligand_name=arguments.ligand_name)

if __name__ == "__main__":
    main()