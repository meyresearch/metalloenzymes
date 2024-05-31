import argparse
import equilibrate
import functions
import shutil
from definitions import KELVIN, FEMTOSECOND, ATM
import BioSimSpace as bss
import os
import meze
import sofra
import solvate


def production(hot_meze, system, outputs, time, timestep=2, temperature=300, pressure=1, restart_write_steps=10000, coordinate_write_steps=10000):

    working_directory = functions.mkdir(outputs + hot_meze.ligand_name + "/")

    if type(timestep) != bss.Types._time.Time:
        timestep = functions.convert_to_units(timestep, FEMTOSECOND)
    if type(temperature) != bss.Types._temperature.Temperature:
        temperature = functions.convert_to_units(temperature, KELVIN)
    if type(pressure) != bss.Types._pressure.Pressure:
        pressure = functions.convert_to_units(pressure, ATM)

    configuration = {}
    restraints_file = None
    if hot_meze.is_metal and hot_meze.restraints:
        configuration = {"nmropt": 1}
        restraints_file = hot_meze.write_restraints_file_0(workdir=working_directory)

    protocol = bss.Protocol.Production(timestep=timestep,
                                       runtime=time,
                                       temperature=temperature,
                                       pressure=pressure,
                                       restart_interval=restart_write_steps,
                                       report_interval=coordinate_write_steps)
    
    checkpoint = hot_meze.equilibration_directory + "/" + hot_meze.ligand_name + "/08_relax/08_relax"

    hot_meze.run(system, protocol, "md", working_directory, configuration, checkpoint, restraints_file)


def main():
    
    parser = argparse.ArgumentParser(description="MEZE: MetalloEnZymE FF-builder for alchemistry")

    parser.add_argument("ligand_name",
                        help="name of ligand",
                        type=str)
    
    parser.add_argument("-n",
                        "--non-metal",
                        dest="non_metal", 
                        help="boolean flag for preparing non-metal proteins",
                        action="store_true")
    
    parser.add_argument("-co",
                        "--cut-off",
                        dest="cut_off",
                        help="cut off distance in Å, used to define zinc coordinating ligands and active site",
                        default=2.6,
                        type=float)
    
    parser.add_argument("-fc0",
                        "--force-constant-0",
                        dest="force_constant_0",
                        help="force constant for model 0",
                        default=100,
                        type=float)    
    
    parser.add_argument("-if",
                        "--input-pdb-file",
                        dest="protein",
                        type=functions.file_exists,
                        help="input pdb file for the metalloenzyme/protein")
    
    parser.add_argument("-g",
                        "--group-name",
                        dest="group_name",
                        help="group name for system, e.g. vim2/kpc2/ndm1/etc.; default taken from input filename",
                        default="meze")

    parser.add_argument("-ppd",
                        "--protein-prep-directory",
                        dest="protein_directory",
                        help="path to protein preparation directory, default: $CWD + /inputs/protein/",
                        default=os.getcwd()+"/inputs/protein/")

    parser.add_argument("-ff",
                        "--force-field",
                        dest="forcefield",
                        help="protein force field, default is ff14SB",
                        default="ff14SB")
    
    parser.add_argument("-w",
                        "--water-model",
                        dest="water_model",
                        help="water model, default is tip3p",
                        default="tip3p")
    
    parser.add_argument("-lpd",
                        "--ligand-prep-directory",
                        dest="ligand_directory",
                        help="path to ligands preparation directory, default: $CWD + /inputs/ligands/",
                        default=os.getcwd()+"/inputs/ligands/")
   
    parser.add_argument("-c",
                        "--ligand-charge",
                        dest="ligand_charge",
                        help="ligand charge",
                        default=0)

    parser.add_argument("-pwd",
                        "--project-working-directory",
                        dest="working_directory",
                        type=functions.path_exists,
                        help="working directory for the project",
                        default=os.getcwd())

    parser.add_argument("-e",
                        "--engine",
                        dest="engine",
                        help="MD engine",
                        default="SOMD")
    
    parser.add_argument("-lf",
                        "--ligand-forcefield",
                        dest="ligand_forcefield",
                        help="ligand force field",
                        default="gaff2")
    
    parser.add_argument("-be",
                        "--box-edges",
                        dest="box_edges",
                        help="box edges in angstroms",
                        type=float,
                        default=20) 
    
    parser.add_argument("-bs",
                        "--box-shape",
                        dest="box_shape",
                        help="box shape",
                        default="cubic") 

    parser.add_argument("-st",
                        "--sampling-time",
                        dest="sampling_time",
                        help="sampling time in nanoseconds",
                        type=float,
                        default=50.0) 

    parser.add_argument("-min",
                        "--minimisation-steps",
                        dest="min_steps",
                        help="number of minimisation steps for equilibration stage",
                        type=int,
                        default=5000)
    
    parser.add_argument("-t",
                        "--equilibration-runtime",
                        dest="equilibration_runtime",
                        help="runtime in ps for equilibration steps",
                        type=float,
                        default=1000) 
    
    parser.add_argument("-wt",
                        "--restraint-weight",
                        dest="restraint_weight",
                        help="positional restraint for equilibrations in kcal/mol*A^-2",
                        type=float,
                        default=100)
    
    parser.add_argument("-ntwr",
                        "--restart-write-steps",
                        dest="restart_write_steps",
                        help="number of steps in which 'restrt' files are written, equivalent to amber ntwr",
                        type=float,
                        default=1000)
    
    parser.add_argument("-ntpr",
                        "--coordinate-write-steps",
                        dest="coordinate_write_steps",
                        help="number of steps in which 'mdcrd' files are written, equivalent to amber ntpr option",
                        type=float,
                        default=1000)
    
    parser.add_argument("-p",
                        "--pressure",
                        dest="pressure",
                        help="simulation pressure, in atm",
                        type=float,
                        default=1)
    
    parser.add_argument("-temp",
                        "--temperature",
                        dest="temperature",
                        help="simulation temperature, in Kelvin",
                        type=float,
                        default=300)
    
    parser.add_argument("-edt",
                        "--equilibration-timestep",
                        dest="equilibration_timestep",
                        help="timestep (in fs) for equilibration simulations",
                        type=float,
                        default=1)    
    
    parser.add_argument("--em-step",
                        dest="emstep",
                        help="Step size for energy minimisation",
                        type=float,
                        default=0.001)
    
    parser.add_argument("--em-tolerance",
                        dest="emtol",
                        help="kJ mol-1 nm-1, Maximum force tolerance for energy minimisation",
                        type=float,
                        default=1000)
    
    parser.add_argument("-sm",
                        "--solvation-method",
                        dest="solvation_method",
                        help="MD engine used for solvation; default is BioSimSpace (GROMACS)",
                        choices=["amber", "gromacs"],
                        default="amber",
                        type=str)
    
    parser.add_argument("-sc",
                        "--solvent-closeness",
                        dest="solvent_closeness",
                        help="control how close, in \AA, solvent atoms can come to solute atoms, default 1.0",
                        type=float, 
                        default=1.0)    
    
    parser.add_argument("-cs",
                        "--cutoff-scheme",
                        dest="cutoff_scheme",
                        help="PME or Reaction Field",
                        default="rf",
                        choices=["pme", "rf"],
                        type=str)
    
    parser.add_argument("--no-restraints",
                        dest="no_restraints",
                        help="do not apply harmonic restraints on metal-coordinating residues",
                        action="store_true")
    
    arguments = parser.parse_args()

    if not arguments.non_metal:
        metal = True
    elif arguments.non_metal:
        metal = False

    if arguments.no_restraints:
        apply_restraints = False
    else:
        apply_restraints = True

    
    if metal:
        network = meze.Meze(protein_file=arguments.protein,
                            cut_off=arguments.cut_off,
                            is_md=True,
                            restraints=apply_restraints,
                            md_input_directory=arguments.working_directory + "md_input_files/",
                            workdir=arguments.working_directory,
                            ligand_path=arguments.ligand_directory,
                            ligand_charge=arguments.ligand_charge,
                            ligand_ff=arguments.ligand_forcefield,
                            group_name=arguments.group_name,        
                            protein_path=arguments.protein_directory,
                            water_model=arguments.water_model,
                            protein_ff=arguments.forcefield,
                            engine=arguments.engine,
                            sampling_time=arguments.sampling_time,
                            box_edges=arguments.box_edges,
                            box_shape=arguments.box_shape,
                            cutoff_scheme=arguments.cutoff_scheme,
                            solvation_method=arguments.solvation_method,
                            solvent_closeness=arguments.solvent_closeness,
                            short_nvt=arguments.equilibration_runtime,
                            nvt=arguments.equilibration_runtime,
                            npt=arguments.equilibration_runtime)
    else:
        network = sofra.Sofra(protein_file=arguments.protein,
                              is_md=True,
                              md_input_directory=arguments.working_directory + "md_input_files/",
                              workdir=arguments.working_directory,
                              ligand_path=arguments.ligand_directory,
                              ligand_charge=arguments.ligand_charge,
                              ligand_ff=arguments.ligand_forcefield,
                              group_name=arguments.group_name,        
                              protein_path=arguments.protein_directory,
                              water_model=arguments.water_model,
                              protein_ff=arguments.forcefield,
                              engine=arguments.engine,
                              sampling_time=arguments.sampling_time,
                              box_edges=arguments.box_edges,
                              box_shape=arguments.box_shape,
                              cutoff_scheme=arguments.cutoff_scheme,
                              solvation_method=arguments.solvation_method,
                              solvent_closeness=arguments.solvent_closeness,
                              short_nvt=arguments.equilibration_runtime,
                              nvt=arguments.equilibration_runtime,
                              npt=arguments.equilibration_runtime)

        
    network.protein_water_complex = network.protein.create_complex()
    network.prepared_protein = network.protein.tleap(network.protein_water_complex)
    network.log_directory = network.create_directory(f"{network.working_directory}/logs/")
    network.protocol_file = network.create_protocol_file()

    # solvate.solvate_bound(network, arguments.ligand_name)

    cold_meze = equilibrate.coldMeze(group_name=network.group_name,
                                     ligand_name=arguments.ligand_name,
                                     prepared=True,
                                     restraints=apply_restraints,
                                     is_metal=metal,
                                     equilibration_directory=network.equilibration_directory,
                                     input_protein_file=network.prepared_protein,
                                     protein_directory=network.protein_path,
                                     ligand_directory=network.ligand_path,
                                     min_steps=network.min_steps,
                                     short_timestep=arguments.equilibration_timestep,
                                     min_dt=network.min_dt,
                                     min_tol=network.min_tol,
                                     temperature=network.temperature,
                                     pressure=network.pressure,
                                     restraint_weight=arguments.restraint_weight,
                                     restart_write_steps=arguments.restart_write_steps,
                                     coordinate_write_steps=arguments.coordinate_write_steps,
                                     afe_input_directory="",
                                     is_md=True,
                                     outputs=network.outputs,
                                     log_directory=network.log_directory,
                                     short_nvt=network.short_nvt, nvt=network.nvt, npt=network.npt)

    hot_meze, relaxed_system = cold_meze.heat_md() 

    production(hot_meze=hot_meze, 
               system=relaxed_system, 
               outputs=hot_meze.outputs,
               time=network.md_time)


if __name__ == "__main__":
    main()
