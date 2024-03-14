import argparse
import equilibrate
import functions
import shutil
from definitions import NANOSECOND, PICOSECOND, KELVIN
import BioSimSpace as bss
import os




def main():
    
    parser = argparse.ArgumentParser(description="MEZE: MetalloEnZymE FF-builder for alchemistry")

    parser.add_argument("ligand_name",
                        help="name of ligand",
                        type=SyntaxWarning)

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
                        default=4) 
    
    parser.add_argument("-r",
                        "--repeats",
                        dest="repeats",
                        help="number of AFE repeat runs",
                        type=int, 
                        default=3)

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

    arguments = parser.parse_args()

    cold_meze = equilibrate.coldMeze(group_name=arguments.group_name,
                                     ligand_name=arguments.ligand_name,
                                     equilibration_directory=arguments.working_directory + "/equilibration/",
                                     input_protein_file=arguments.protein,
                                     protein_directory=arguments.protein_directory,
                                     ligand_directory=arguments.ligand_directory,
                                     min_steps=arguments.min_steps,
                                     short_nvt=arguments.equilibration_runtime,
                                     short_timestep=arguments.equilibration_timestep,
                                     nvt=arguments.equilibration_runtime,
                                     npt=arguments.equilibration_runtime,
                                     min_dt=arguments.emstep,
                                     min_tol=arguments.emtol,
                                     temperature=arguments.temperature,
                                     pressure=arguments.pressure,
                                     restraint_weight=arguments.restraint_weight,
                                     restart_write_steps=arguments.restart_write_steps,
                                     coordinate_write_steps=arguments.coordinate_write_steps)

    cold_meze.heat_md() 


if __name__ == "__main__":
    main()


# To be depracated:
def equilibration(system, ligand_name, minimised_system, nonbonded_cut_off=8.0, dt=0.001, runtime=1, start_temperature=100): 
        
    directory = system.equilibration_directory+f"{ligand_name}/"
    directories = lambda step: functions.mkdir(directory + step)
    
    heat_02_dir = directories("02_heat")
    relax_03_dir = directories("03_relax")
    lower_04_dir = directories("04_lower")
    bb_min_05_dir = directories("05_bb_min")
    relax_06_dir = directories("06_relax")
    reduce_07_dir = directories("07_reduce")
    continue_08_dir = directories("08_continue")
    relax_09_dir = directories("09_relax")

    template_restraints_file = system.write_restraints_file_0()

    restraints_files = []
    for directory in [heat_02_dir, relax_03_dir, lower_04_dir, bb_min_05_dir, relax_06_dir, reduce_07_dir, continue_08_dir, relax_09_dir]:
        restraints_files.append(shutil.copy(template_restraints_file, directory).split("/")[-1])

    max_cycles = system.min_steps        
    output_frequency = max_cycles // 10

    heat_02_options = {"ioutfm": 1,
                        "cut": nonbonded_cut_off,
                        "tol": 0.00001,
                        "nscm": 0,
                        "barostat": 2,
                        "taup": 1.0,
                        "gamma_ln": 1.0,
                        "ntp": 0,
                        "ntb": 1, 
                        "ntc": 2,
                        "ntf": 2,   
                        "nmropt": 1}
    
    runtime_ns = functions.convert_to_units(runtime, NANOSECOND)
    namelist = ["&wt TYPE='DUMPFREQ', istep1=1 /"]
    heat_02_protocol = bss.Protocol.Equilibration(timestep=dt*PICOSECOND, 
                                                runtime=runtime_ns, 
                                                temperature_start=start_temperature*KELVIN, 
                                                temperature_end=system.temperature, 
                                                report_interval=output_frequency, 
                                                restart_interval=output_frequency,
                                                restraint="heavy",
                                                force_constant=system.force_constant_0)
    amber_home = os.environ["AMBERHOME"]
    heat_02_process = bss.Process.Amber(system=minimised_system, 
                                        protocol=heat_02_protocol, 
                                        name="02_heat", 
                                        work_dir=heat_02_dir, 
                                        extra_options=heat_02_options,
                                        extra_lines=namelist,
                                        exe=amber_home + "/bin/pmemd.cuda")
    heat_02_config = heat_02_dir + "/*.cfg"
    heat_02_config_file = functions.get_files(heat_02_config)[0]        

    with open(heat_02_config_file, "a") as file:
        file.write("\n")
        file.write(f"DISANG={restraints_files[0]}\n")
        file.write(f"DUMPAVE=distances.out\n")

    heat_02_process.start()
    heat_02_process.wait()
    if heat_02_process.isError():
        print(heat_02_process.stdout())
        print(heat_02_process.stderr())
        raise bss._Exceptions.ThirdPartyError("The process exited with an error!")

    heat_02_system = heat_02_process.getSystem()

    relax_03_options = {"ioutfm": 1,
                        "cut": nonbonded_cut_off,
                        "nscm": 0,
                        "barostat": 2,
                        "ntp": 1,
                        "taup": 1.0,
                        "gamma_ln": 1.0,
                        "ntp": 1,
                        "ntb": 2, # constant pressure 
                        "irest": 1,
                        "ntx": 5,
                        "ntr": 1,
                        "nmropt": 1}
    
    relax_03_protocol = bss.Protocol.Equilibration(timestep=dt*PICOSECOND, 
                                                runtime=runtime_ns, 
                                                temperature_start=0.0*KELVIN, 
                                                temperature_end=system.temperature, 
                                                report_interval=output_frequency, 
                                                restart_interval=output_frequency,
                                                restraint="heavy",
                                                force_constant=system.force_constant_0)
    
    relax_03_process = bss.Process.Amber(system=heat_02_system,  
                                         protocol=relax_03_protocol, 
                                         name="03_relax", 
                                         work_dir=relax_03_dir, 
                                         extra_options=relax_03_options,
                                         extra_lines=namelist,
                                         exe=amber_home + "/bin/pmemd.cuda")
    

    relax_03_config = relax_03_dir + "/*.cfg"
    relax_03_config_file = functions.get_files(relax_03_config)[0]

    with open(relax_03_config_file, "r") as file:
        config = file.readlines()

    new_config = [line for line in config if "tempi=" not in line and "TEMP0" not in line]
    
    with open(relax_03_config_file, "w") as file:
        file.writelines(new_config)
        file.write("\n")
        file.write(f"DISANG={restraints_files[1]}\n")
        file.write(f"DUMPAVE=distances.out\n")

    relax_03_process.start()
    relax_03_process.wait()
    if relax_03_process.isError():
        print(relax_03_process.stdout())
        print(relax_03_process.stderr())
        raise bss._Exceptions.ThirdPartyError("The process exited with an error!")

    relax_03_system = relax_03_process.getSystem()

    lower_04_options = relax_03_options
    lower_04_protocol = bss.Protocol.Equilibration(timestep=dt*PICOSECOND, 
                                                    runtime=runtime_ns, 
                                                    temperature_start=0.0*KELVIN, 
                                                    temperature_end=system.temperature, 
                                                    report_interval=output_frequency, 
                                                    restart_interval=output_frequency,
                                                    restraint="heavy",
                                                    force_constant=system.force_constant_0 / 10)

    lower_04_process = bss.Process.Amber(system=relax_03_system,  
                                         protocol=lower_04_protocol, 
                                         name="04_lower", 
                                         work_dir=lower_04_dir, 
                                         extra_options=lower_04_options,
                                         extra_lines=namelist,
                                         exe=amber_home + "/bin/pmemd.cuda")
    
    lower_04_config = lower_04_dir + "/*.cfg"
    lower_04_config_file = functions.get_files(lower_04_config)[0]

    with open(lower_04_config_file, "r") as file:
        config = file.readlines()

    new_config = [line for line in config if "tempi=" not in line and "TEMP0" not in line]
    
    with open(lower_04_config_file, "w") as file:
        file.writelines(new_config)
        file.write("\n")
        file.write(f"DISANG={restraints_files[2]}\n")
        file.write(f"DUMPAVE=distances.out\n")

    lower_04_process.start()
    lower_04_process.wait()
    if lower_04_process.isError():
        print(lower_04_process.stdout())
        print(lower_04_process.stderr())
        raise bss._Exceptions.ThirdPartyError("The process exited with an error!")

    lower_04_system = lower_04_process.getSystem()


    max_cycles = system.min_steps 
    n_steep_cycles = max_cycles // 33

    bb_min_05_options = {"ntmin": 1,
                            "ntpr": output_frequency,
                            "ntwx": output_frequency,
                            "ntpr": output_frequency,
                            "ntp": 0,
                            "ioutfm": 1,
                            "cut": nonbonded_cut_off, 
                            "ncyc": n_steep_cycles,
                            "ntc": 2, # constrain H bonds
                            "ntf": 2, # exlcude H bonds from force calc.
                            "ntb": 1, # constant volume
                            "ntr": 1,
                            "nmropt": 1}   
    
    bb_min_05_protocol = bss.Protocol.Minimisation(steps=max_cycles,
                                                    restraint="backbone",
                                                    force_constant=system.force_constant_0 / 10)

    bb_min_05_process = bss.Process.Amber(system=lower_04_system, 
                                                protocol=bb_min_05_protocol, 
                                                name="05_bb_min", 
                                                work_dir=bb_min_05_dir, 
                                                extra_options=bb_min_05_options,
                                                extra_lines=namelist,
                                                exe=amber_home + "/bin/pmemd.cuda")
    
    bb_min_05_config = bb_min_05_dir + "/*.cfg"
    bb_min_05_config_file = functions.get_files(bb_min_05_config)[0]

    with open(bb_min_05_config_file, "r") as file:
        config = file.readlines()
    new_config = [line for line in config if "ig=" not in line]
    
    with open(bb_min_05_config_file, "w") as file:
        file.writelines(new_config)
        file.write("\n")
        file.write(f"DISANG={restraints_files[3]}\n")
        file.write(f"DUMPAVE=distances.out\n")

    bb_min_05_process.start()
    bb_min_05_process.wait()
    if bb_min_05_process.isError():
        print(bb_min_05_process.stdout())
        print(bb_min_05_process.stderr())
        raise bss._Exceptions.ThirdPartyError("The process exited with an error!")
    
    bb_min_05_system = bb_min_05_process.getSystem()        

    relax_06_options = lower_04_options
    relax_06_protocol = bss.Protocol.Equilibration(timestep=dt*PICOSECOND, 
                                                runtime=runtime_ns, 
                                                temperature_start=system.temperature, 
                                                temperature_end=system.temperature, 
                                                report_interval=output_frequency, 
                                                restart_interval=output_frequency,
                                                restraint="backbone",
                                                force_constant=system.force_constant_0 / 10)

    relax_06_process = bss.Process.Amber(system=bb_min_05_system,  #CHANGE
                                        protocol=relax_06_protocol, 
                                        name="06_relax", 
                                        work_dir=relax_06_dir, 
                                        extra_options=relax_06_options,
                                        extra_lines=namelist,
                                        exe=amber_home + "/bin/pmemd.cuda")        

    relax_06_config = relax_06_dir + "/*.cfg"
    relax_06_config_file = functions.get_files(relax_06_config)[0]   

    with open(relax_06_config_file, "r") as file:
        config = file.readlines()

    new_config = [line for line in config if "tempi=" not in line and "TEMP0" not in line]
            

    with open(relax_06_config_file, "w") as file:
        file.writelines(new_config)
        file.write("\n")
        file.write(f"DISANG={restraints_files[4]}\n")
        file.write(f"DUMPAVE=distances.out\n")

    relax_06_process.start()
    relax_06_process.wait()
    if relax_06_process.isError():
        print(relax_06_process.stdout())
        print(relax_06_process.stderr())
        raise bss._Exceptions.ThirdPartyError("The process exited with an error!")        

    relax_06_system = relax_06_process.getSystem()

    reduce_07_options = relax_06_options

    reduce_07_protocol = bss.Protocol.Equilibration(timestep=dt*PICOSECOND, 
                                                runtime=runtime_ns, 
                                                temperature_start=system.temperature, 
                                                temperature_end=system.temperature, 
                                                report_interval=output_frequency, 
                                                restart_interval=output_frequency,
                                                restraint="backbone",
                                                force_constant=system.force_constant_0 / system.force_constant_0)
    
    reduce_07_process = bss.Process.Amber(system=relax_06_system,  #CHANGE
                                        protocol=reduce_07_protocol, 
                                        name="07_reduce", 
                                        work_dir=reduce_07_dir, 
                                        extra_options=reduce_07_options,
                                        extra_lines=namelist,
                                        exe=amber_home + "/bin/pmemd.cuda")

    reduce_07_config = reduce_07_dir + "/*.cfg"
    reduce_07_config_file = functions.get_files(reduce_07_config)[0]          

    with open(reduce_07_config_file, "r") as file:
        config = file.readlines()

    new_config = [line for line in config if "tempi=" not in line and "TEMP0" not in line]
            
    with open(reduce_07_config_file, "w") as file:
        file.writelines(new_config)
        file.write("\n")
        file.write(f"DISANG={restraints_files[5]}\n")
        file.write(f"DUMPAVE=distances.out\n")

    reduce_07_process.start()
    reduce_07_process.wait()
    if relax_06_process.isError():
        print(reduce_07_process.stdout())
        print(reduce_07_process.stderr())
        raise bss._Exceptions.ThirdPartyError("The process exited with an error!")        

    reduce_07_system = reduce_07_process.getSystem()

    continue_08_options = reduce_07_options

    continue_08_protocol = bss.Protocol.Equilibration(timestep=dt*PICOSECOND, 
                                                      runtime=runtime_ns, 
                                                      temperature_start=system.temperature, 
                                                      temperature_end=system.temperature, 
                                                      report_interval=output_frequency, 
                                                      restart_interval=output_frequency,
                                                      restraint="backbone",
                                                      force_constant=(system.force_constant_0 / system.force_constant_0) / 10)
    
    continue_08_process = bss.Process.Amber(system=reduce_07_system,
                                            protocol=continue_08_protocol,
                                            name="08_continue",
                                            work_dir=continue_08_dir,
                                            extra_options=continue_08_options,
                                            extra_lines=namelist, 
                                            exe=amber_home + "/bin/pmemd.cuda")
    
    continue_08_config = continue_08_dir + "/*.cfg"
    continue_08_config_file = functions.get_files(continue_08_config)

    with open(continue_08_config_file, "r") as file:
        config = file.readlines()

    new_config = [line for line in config if "tempi=" not in line and "TEMP0" not in line]
                
    with open(continue_08_config_file, "w") as file:
        file.writelines(new_config)
        file.write("\n")
        file.write(f"DISANG={restraints_files[6]}\n")
        file.write(f"DUMPAVE=distances.out\n")

    continue_08_process.start()
    continue_08_process.wait()
    if continue_08_process.isError():
        print(continue_08_process.stdout())
        print(continue_08_process.stderr())
        raise bss._Exceptions.ThirdPartyError("The process exited with an error!")        

    continue_08_system = continue_08_process.getSystem()        

    relax_09_options = {"ioutfm": 1,
                        "cut": nonbonded_cut_off,
                        "tol": 0.00001,
                        "nscm": max_cycles,
                        "barostat": 2,
                        "ntp": 1,
                        "taup": 1.0,
                        "gamma_ln": 1.0,
                        "ntb": 2, 
                        "ntc": 2,
                        "ntf": 2,
                        "iwrap": 0,  
                        "irest": 1,
                        "ntx": 5, 
                        "nmropt": 1}

    relax_09_protocol = bss.Protocol.Equilibration(timestep=dt*PICOSECOND, 
                                                runtime=runtime_ns, 
                                                temperature=system.temperature, 
                                                report_interval=output_frequency, 
                                                restart_interval=output_frequency)

    relax_09_process = bss.Process.Amber(system=continue_08_system,
                                            protocol=relax_09_protocol,
                                            name="09_relax",
                                            work_dir=relax_09_dir,
                                            extra_options=relax_09_options,
                                            extra_lines=namelist,
                                            exe=amber_home + "/bin/pmemd.cuda")

    relax_09_config = relax_09_dir + "/*.cfg"
    relax_09_config_file = functions.get_files(relax_09_config)[0]        

    with open(relax_09_config_file, "a") as file:
        file.write("\n")
        file.write(f"DISANG={restraints_files[7]}\n")
        file.write(f"DUMPAVE=distances.out\n")        

    relax_09_process.start()
    relax_09_process.wait()
    if relax_09_process.isError():
        print(relax_09_process.stdout())
        print(relax_09_process.stderr())
        raise bss._Exceptions.ThirdPartyError("The process exited with an error!")
    
    equilibrated_system = relax_09_process.getSystem()

    return equilibrated_system  


def minimisation(system, ligand_name, nonbonded_cut_off=10.0):

    directory = functions.mkdir(system.equilibration_directory+f"{ligand_name}/")
    files = functions.get_files(f"{system.protein_path}/bound_{ligand_name}_solvated.*")
    solvated_system = bss.IO.readMolecules(files)
    directories = lambda step: functions.mkdir(directory + step)
    min_dir = directories("01_min/")

    restraints_file = system.write_restraints_file_0(workdir=min_dir)

    max_cycles = system.min_steps        
    n_steep_cycles = max_cycles // 20
    output_frequency = max_cycles // 20

    minimisation_options = {"ntmin": 1,
                            "ntpr": output_frequency,
                            "ntwx": output_frequency,
                            "ntpr": output_frequency,
                            "ntp": 0,
                            "ioutfm": 1,
                            "cut": nonbonded_cut_off, 
                            "ncyc": n_steep_cycles,
                            "ntc": 2, # constrain H bonds
                            "ntf": 2, # exlcude H bonds from force calc.
                            "ntb": 1, # constant volume
                            "nmropt": 1}       
    
    namelist = ["&wt TYPE='DUMPFREQ', istep1=1 /"]

    minimisation_protocol = bss.Protocol.Minimisation(steps=max_cycles,
                                                      restraint="heavy",
                                                      force_constant=system.force_constant_0)

    amber_home = os.environ["AMBERHOME"]
    minimisation_process = bss.Process.Amber(system=solvated_system, 
                                             protocol=minimisation_protocol, 
                                             name="01_min", 
                                             work_dir=min_dir, 
                                             extra_options=minimisation_options,
                                             extra_lines=namelist,
                                             exe=amber_home + "/bin/pmemd.cuda")
    
    min_config = min_dir + "/*.cfg"
    config_file = functions.get_files(min_config)[0]

    with open(config_file, "r") as file:
        config = file.readlines()
    new_config = [line for line in config if "ig=" not in line]
    
    with open(config_file, "w") as file:
        file.writelines(new_config)
        file.write("\n")
        file.write(f"DISANG={functions.get_filename(restraints_file)}.RST\n")
        file.write(f"DUMPAVE=distances.out\n")

    minimisation_process.start()
    minimisation_process.wait()
    if minimisation_process.isError():
        print(minimisation_process.stdout())
        print(minimisation_process.stderr())
        raise bss._Exceptions.ThirdPartyError("The process exited with an error!")
    system = minimisation_process.getSystem()
    return system

def production(system, ligand_name, equilibrated_system, nonbonded_cut_off=9.0, dt=0.002, runtime=50):

    directory = system.output_directory+f"{ligand_name}/"
    production_directory = functions.mkdir(directory)

    template_restraints_file = system.write_restraints_file_0()
    restraints_file = shutil.copy(template_restraints_file, production_directory).split("/")[-1]
    output_frequency = runtime // dt
    production_options = {"ioutfm": 1,
                            "cut": nonbonded_cut_off,
                            "gamma_ln": 1.0,
                            "ntp": 1,
                            "ntb": 2, 
                            "ntc": 2,
                            "ntf": 2, 
                            "iwrap": 1,  
                            "nmropt": 1}
    
    runtime_ns = functions.convert_to_units(runtime, NANOSECOND)
    namelist = ["&wt TYPE='DUMPFREQ', istep1=1 /"]

    protocol = bss.Protocol.Production(timestep=dt*PICOSECOND,
                                        runtime=runtime_ns,
                                        temperature=system.temperature,
                                        report_interval=output_frequency,
                                        restart_interval=output_frequency,
                                        restart=True)
    
    process = bss.Process.Amber(system=equilibrated_system, 
                                protocol=protocol, 
                                name="md",
                                work_dir=production_directory,
                                extra_options=production_options,
                                extra_lines=namelist)
    config = production_directory + "/*.cfg"
    config_file = functions.get_files(config)[0]

    with open(config_file, "a") as file:
        file.write("\n")
        file.write(f"DISANG={restraints_file}\n")
        file.write(f"DUMPAVE=distances.out\n")

    process.start()
    process.wait()
    if process.isError():
        print(process.stdout())
        print(process.stderr())
        raise bss._Exceptions.ThirdPartyError("The process exited with an error!")

def run_md(system, ligand_name):

    minimised_system = minimisation(system, ligand_name)
    equilibrated_system = equilibration(system, ligand_name, minimised_system)
    production(system, ligand_name, equilibrated_system)
