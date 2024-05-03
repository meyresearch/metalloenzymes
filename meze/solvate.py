import argparse
import BioSimSpace as bss
import Ligand
from definitions import ANGSTROM
import sofra 
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
    Ligand:
        (solvated) Ligand object whose file attribute is the prm7 and rst7 files 
    """ 
    method = network.solvation_method
    ligand = network.get_ligand_by_name(ligand_name)
    ligand_savename = network.ligand_path + ligand_name + "_solvated"
    print(f"Solvating unbound ligand {ligand_name}")

    if method == "gromacs":
        ligand_parameters = ligand.parameterise(network.ligand_forcefield, network.ligand_charge)
        unbound_box, unbound_box_angles = network.create_box(ligand_parameters)
        # bss.IO.saveMolecules(f"{network.ligand_path}/test", ligand_parameters, "pdb")
        solvated_molecule = bss.Solvent.solvate(model=network.protein.water_model, 
                                                molecule=ligand_parameters, 
                                                box=unbound_box,
                                                angles=unbound_box_angles,
                                                work_dir=network.ligand_path)
        solvated_files = bss.IO.saveMolecules(ligand_savename, solvated_molecule, ["PRM7", "RST7"])

    elif method == "amber":
        
        solvate_directory = network.create_directory(f"{network.ligand_path}/solvate_{ligand_name}/")

        if not os.path.isfile(f"{solvate_directory}/{ligand_name}.mol2"):
            parameterised_ligand = ligand.antechamber(charge=network.ligand_charge, path=solvate_directory, atom_type=network.ligand_forcefield)
            ligand_mol2_file = parameterised_ligand.file
        else:
            ligand_mol2_file = functions.get_files(f"{solvate_directory}/{ligand_name}.mol2")[0]
        
        if not os.path.isfile(f"{solvate_directory}/{ligand_name}.frcmod"):
            ligand_frcmod_file = ligand.parmcheck(path=solvate_directory)
        else:
            ligand_frcmod_file = functions.get_files(f"{solvate_directory}/{ligand_name}.frcmod")[0]
        
        ligand_pdb = f"{solvate_directory}/{ligand_name}.pdb"
        ligand_tleap_input_file = solvate_directory + f"/tleap_genpdb_{ligand_name}.in"
        ligand_tleap_output_file = solvate_directory + f"/tleap_genpdb_{ligand_name}.out"
        
        if not os.path.isfile(ligand_pdb):
            tleap_pdb_file = f"{ligand_name}_tleap.pdb"
            with open(ligand_tleap_input_file, "w") as tleap_in:
                tleap_in.write(f"source leaprc.{network.ligand_forcefield}\n")
                tleap_in.write(f"source leaprc.water.{network.water_model}\n")
                tleap_in.write(f"lig = loadmol2 {ligand_mol2_file}\n")
                tleap_in.write(f"loadamberparams {ligand_frcmod_file}\n")
                tleap_in.write(f"savepdb lig {tleap_pdb_file}\n")
                tleap_in.write("quit")
            work_dir = os.getcwd()
            os.chdir(solvate_directory)
            os.system(f"tleap -s -f {ligand_tleap_input_file} > {ligand_tleap_output_file}")
            os.system(f"pdb4amber -i {tleap_pdb_file} -o {ligand_pdb}")
            os.chdir(work_dir)

        tleap_input_file = solvate_directory + f"/tleap_{ligand_name}_solvate.in"
        tleap_output_file = solvate_directory + f"/tleap_{ligand_name}_solvate.out"

        if not os.path.isfile(f"{ligand_savename}.prm7") or not os.path.isfile(f"{ligand_savename}.rst7"):
            with open(tleap_input_file, "w") as tleap_in:
                tleap_in.write(f"source leaprc.{network.ligand_forcefield}\n")
                tleap_in.write(f"source leaprc.water.{network.water_model}\n")
                if network.water_model == "tip3p":
                    tleap_in.write(f"loadamberparams frcmod.ions1lm_126_tip3p\n")
                tleap_in.write("\n")
                tleap_in.write(f"lig = loadmol2 {ligand_mol2_file}\n")
                tleap_in.write(f"loadamberparams {ligand_frcmod_file}\n")
                tleap_in.write("\n")
                box_shape_three_letter = network.box_shape[:3]
                tleap_in.write(f"solvate{box_shape_three_letter} lig {network.water_model.upper()}BOX {network.box_edges} {network.solvent_closeness}\n")
                tleap_in.write(f"addions lig Na+ 0\n")
                tleap_in.write(f"addions lig Cl- 0\n")
                tleap_in.write("\n")
                tleap_in.write(f"saveamberparm lig {ligand_savename}.prm7 {ligand_savename}.rst7\n")
                tleap_in.write(f"quit\n")

            work_dir = os.getcwd()
            os.chdir(solvate_directory)
            os.system(f"tleap -s -f {tleap_input_file} > {tleap_output_file}")
            os.chdir(work_dir)
            
        solvated_files = functions.get_files(ligand_savename + ".*")
    
    return Ligand.Ligand(file=solvated_files, parameterised=True)



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
    method = network.solvation_method
    ligand = network.get_ligand_by_name(ligand_name)
    system_savename = network.protein_path + "bound_" + ligand_name + "_solvated"
    print(f"Solvating bound ligand {ligand_name}")

    if method == "gromacs":
        ligand_parameters = ligand.parameterise(network.ligand_forcefield, network.ligand_charge)    
        protein = network.protein.get_molecule()
        system_parameters = ligand_parameters + protein
        bound_box, bound_box_angles = network.create_box(system_parameters)
        solvated_molecules = bss.Solvent.solvate(model=network.protein.water_model,
                                                    molecule=system_parameters,
                                                    box=bound_box,
                                                    angles=bound_box_angles)
        
        
        solvated_files = bss.IO.saveMolecules(system_savename, solvated_molecules, ["PRM7", "RST7"])

    elif method == "amber":
        
        solvate_directory = network.create_directory(f"{network.protein_path}/solvate_{ligand_name}_bound/")
        
        if not os.path.isfile(f"{network.ligand_path}/solvate_{ligand_name}/{ligand_name}.mol2"):
            parameterised_ligand = ligand.antechamber(charge=network.ligand_charge, path=solvate_directory, atom_type=network.ligand_forcefield)
            ligand_mol2_file = parameterised_ligand.file
        else:
            ligand_mol2_file = functions.get_files(f"{network.ligand_path}/solvate_{ligand_name}/{ligand_name}.mol2")[0]
        
        if not os.path.isfile(f"{network.ligand_path}/solvate_{ligand_name}/{ligand_name}.frcmod"):
            ligand_frcmod_file = ligand.parmcheck(path=solvate_directory)
        else:
            ligand_frcmod_file = functions.get_files(f"{network.ligand_path}/solvate_{ligand_name}/{ligand_name}.frcmod")[0]

        tleap_input_file = solvate_directory + f"/tleap_bound_{ligand_name}_solvate.in"
        tleap_output_file = solvate_directory + f"/tleap_bound_{ligand_name}_solvate.out"   

        protein_pdb = functions.get_files(functions.file_exists(f"{network.protein_path}/{network.group_name}.pdb"))[0]

        if not os.path.isfile(f"{system_savename}.prm7") or not os.path.isfile(f"{system_savename}.rst7"):
            with open(tleap_input_file, "w") as tleap_in:
                tleap_in.write(f"source oldff/leaprc.{network.protein_forcefield}\n")
                tleap_in.write(f"source leaprc.{network.ligand_forcefield}\n")
                tleap_in.write(f"source leaprc.water.{network.water_model}\n")
        
                if network.water_model == "tip3p":
                    tleap_in.write(f"loadamberparams frcmod.ions1lm_126_tip3p\n")
                tleap_in.write(f"loadamberparams {ligand_frcmod_file}\n")   
                
                tleap_in.write("\n")
                tleap_in.write(f"lig = loadmol2 {ligand_mol2_file}\n")
                tleap_in.write("\n")
                tleap_in.write(f"protein = loadpdb {protein_pdb}\n")

                tleap_in.write("complex = combine {protein lig}\n")
                tleap_in.write(f"savepdb complex {ligand_name}_complex_dry.pdb\n")
                tleap_in.write("check complex\n")

                box_shape_three_letter = network.box_shape[:3]
                tleap_in.write(f"solvate{box_shape_three_letter} complex {network.water_model.upper()}BOX {network.box_edges} {network.solvent_closeness}\n")
                tleap_in.write(f"addions2 complex Na+ 0\n")
                tleap_in.write(f"addions2 complex Cl- 0\n")
                tleap_in.write("\n")
                tleap_in.write(f"savepdb complex {ligand_name}_complex_solvated.pdb\n")
                tleap_in.write(f"saveamberparm complex {system_savename}.prm7 {system_savename}.rst7\n")
                tleap_in.write(f"quit\n")
            
            work_dir = os.getcwd()
            os.chdir(solvate_directory)
            os.system(f"tleap -s -f {tleap_input_file} > {tleap_output_file}")
            os.chdir(work_dir)

        solvated_files = functions.get_files(system_savename + ".*")
            
    return Ligand.Ligand(file=solvated_files, parameterised=True)


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
                        help="protocol file",
                        type=str,
                        default=os.getcwd() + "/afe/protocol.dat")
    
    arguments = parser.parse_args()
    
    protocol = functions.input_to_dict(arguments.protocol_file)

    network = sofra.Sofra(prepared=True,
                          workdir=protocol["project directory"],
                          afe_input_path=protocol["afe input directory"],
                          equilibration_path=protocol["equilibration directory"],
                          outputs=protocol["outputs"],
                          log_directory=protocol["log directory"],
                          ligand_path=protocol["ligand directory"],
                          group_name=protocol["group name"],
                          protein_file=protocol["prepared protein file"],
                          protein_path=protocol["protein directory"],
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
                          pressure=protocol["pressure"],
                          solvation_method=protocol["solvation method"],
                          solvent_closeness=protocol["solvent closeness"])

    solvate_unbound(network=network, ligand_name=arguments.ligand_name)
    solvate_bound(network=network, ligand_name=arguments.ligand_name)

if __name__ == "__main__":
    main()