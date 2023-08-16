"""
Alchemical Free Energy class for free energy simulations
"""
import pathlib
import csv
import BioSimSpace as bss
from definitions import ANGSTROM, PICOSECOND, NANOSECOND, KELVIN, ATM
import functions 


def run_process(system, protocol, process, working_directory, configuration=None):
    """
    Run a Gromacs minimisation or equilibration process 
    Adapted from https://tinyurl.com/BSSligprep

    Parameters:
    -----------
    system: bss.System
        run system
    protocol: bss.Protocol 
        minimisation or equilibration
    process: name 
        process name for saving process output
    working_directory: str
        save output into this directory

    Return:
    -------
    system: bss.System
        equilibrated or minimised system
    """
    process = bss.Process.Gromacs(system, protocol, name=process, work_dir=working_directory)
    config = process.getConfig()
    if configuration:
        for setting in configuration:
            key = setting.split()[0]
            try:
                index = [i for i, string in enumerate(config) if key in string][0]
                config[index] = setting
                process.setConfig(config)
            except IndexError:
                process.addToConfig(setting)
                config = process.getConfig()
    process.setArg("-ntmpi", 1)
    process.start()
    process.wait()
    if process.isError():
        print(process.stdout())
        print(process.stderr())
        raise bss._Exceptions.ThirdPartyError("The process exited with an error!")
    system = process.getSystem()
    return system


class Simulation(object):
    """
    _summary_

    Attributes:
    -----------
    object: 
        _description_

    Methods:
    -------
    """
    def __init__(self, path, engine="SOMD", sampling_time=4, box_edges=20, box_shape="cubic"):
        """
        Class constructor
        """
        self.path = functions.path_exists(path)
        self.engine = engine
        self.time = functions.convert_to_units(sampling_time, NANOSECOND)
        self.box_shape = box_shape
        self.box_edges = box_edges
        self.afe_dir = self.create_directory(self.path + "/afe/")

    def create_directory(self, directory):
        """
        Create AFE working directory in path.

        Parameters:
        -----------
        name: str
            name of new directory
        Return:
        -------
        directory: str
            full path to afe directory
        """
        try:
            pathlib.Path(directory).mkdir(parents=False, exist_ok=False)
        except FileNotFoundError as e:
            print(f"Could not create directory {directory}. Pathlib raised error: {e}")
        except FileExistsError as e:
            print(f"Could not create directory {directory}. Pathlib raised error: {e}")
        return directory        
    
    
    def create_dat_file(self, Protein, Network):
        """
        Create protocol.dat file for AFE runs

        Parameters:
        -----------
        Protein: Protein
            Protein class object
        Ligands: Ligands
            Ligands class object

        Return:
        -------
        protocol_file: str
            protocol datafile
        """
        protocol = [f"ligand forcefield = {Network.forcefield}", 
                    f"protein forcefield = {Protein.forcefield}", 
                    f"solvent = {Protein.water_model}", 
                    f"box edges = {self.box_edges}*angstrom", 
                    f"box shape = {self.box_shape}", 
                    f"protocol = default",
                    f"sampling = {self.time}*ns",
                    f"engine = {self.engine}"]
        self.protocol_file = self.home + "/protocol.dat"

        with open(self.protocol_file, "w") as file:
            writer = csv.writer(file)
            for protocol_line in protocol:
                writer.writerow([protocol_line])
        return self.protocol_file
    

    def create_box(self, molecule):
        """
        Create a bss.Box object for solvation.

        Parameters:
        -----------
        molecule: 
            bss.Molecule: usually either a protein or a ligand

        Return:
        -------
        tuple: 
            bss.Box and angles
        """

        box_min, box_max = molecule.getAxisAlignedBoundingBox()
        box_size = [y - x for x, y in zip(box_min, box_max)]
        box_area = [x + int(self.box_edges) * ANGSTROM for x in box_size]
        self.box, self.box_angles = None, None
        if self.box_shape == "cubic":
            self.box, self.box_angles = bss.Box.cubic(max(box_area))
        elif self.box_shape == "rhombicDodecahedronHexagon":
            self.box, self.box_angles = bss.Box.rhombicDodecahedronHexagon(max(box_area))
        elif self.box_shape == "rhombicDodecahedronSquare":
            self.box, self.box_angles = bss.Box.rhombicDodecahedronSquare(max(box_area))
        elif self.box_shape == "truncatedOctahedron":
            self.box, self.box_angles = bss.Box.truncatedOctahedron(max(box_area))
        else:
            print(f"Box shape {self.box_shape} not supported.")
        return self.box, self.box_angles


class EquilibrationSimulation(Simulation):
    def __init__(self, path, engine="SOMD", sampling_time=4, box_edges=20, box_shape="cubic"):
        super().__init__(path, engine, sampling_time, box_edges, box_shape)
        
        self.minimisation_steps = 500
        self.minimisation_stepsize = 0.01
        self.minimisation_tolerance = 1000
        # self.short_nvt_runtime = functions.convert_to_units(5, PICOSECOND)
        # self.nvt_runtime = functions.convert_to_units(50, PICOSECOND)
        # self.npt_runtime = functions.convert_to_units(200, PICOSECOND)
        # self.equilibration_dir = self.create_directory(self.path + "/equilibration/")
        self.start_temperature = functions.convert_to_units(0, KELVIN)
        self.end_temperature = functions.convert_to_units(300, KELVIN)
        self.temperature = None
        self.pressure = None
        self.restraints = None
        self.configuration = None
        self.name = ""


    def minimise(self, system):
        """
        Minimise the system using Gromacs

        Parameters:
        -----------
        system: bss.System
            system to be minimised
        working_directory: str
            current working dir

        Return:
        -------
        minimsed_system: bss.System
            minimised system
        """
        print("Minimisation")
        working_directory = self.create_directory(self.path + "/min/")
        protocol = bss.Protocol.Minimisation(steps=self.minimisation_steps)
        self.configuration = [f"emstep = {self.minimisation_stepsize}", f"emtol = {self.minimisation_tolerance}"]
        minimised_system = run_process(system, protocol, "min", working_directory)
        return minimised_system


    def equilibrate(self, system):
        """
        Run NVT or NPT equilibration

        Parameters:
        -----------
        system: bss.System
            system to be equilibrated

        Return:
        -------
        equilibrated_system: bss.System
            equilibrated system
        """
        working_directory = self.create_directory(self.path + "/" + self.name)
        print(f"{self.name.upper()}")

        # set these in the constructor when defining object
        # if temperature:
        #     temperature = functions.convert_to_units(temperature, KELVIN)
        # if pressure:
        #     pressure = functions.convert_to_units(pressure, ATM)

        protocol = bss.Protocol.Equilibration(runtime=self.time,
                                              temperature_start=self.start_temperature,
                                              temperature_end=self.end_temperature,
                                              temperature=self.temperature,
                                              pressure=self.pressure,
                                              restraint=self.restraints)
        equilibrated_system = run_process(system, protocol, self.name, working_directory, self.configuration)
        return equilibrated_system


    def set_minimisation_steps(self, min_steps):
        self.minimisation_steps = min_steps
    def set_minimisation_stepsize(self, min_stepsize):
        self.minimisation_stepsize = min_stepsize
    def set_minimisation_tolerance(self, min_tolerance):
        self.minimisation_tolerance = min_tolerance
    def set_runtime(self, runtime):
        self.time = functions.convert_to_units(runtime, PICOSECOND)
    # def set_nvt(self, nvt):
    #     self.nvt_runtime = functions.convert_to_units(nvt, PICOSECOND)
    # def set_npt(self, npt):
    #     self.npt_runtime = functions.convert_to_units(npt, PICOSECOND)
    # def set_equilibration_dir(self, directory):
    #     self.equilibration_dir = self.create_directory(directory)
    def set_start_temperature(self, temperature):
        self.start_temperature = functions.convert_to_units(temperature, KELVIN)
    def set_end_temperature(self, temperature):
        self.end_temperature = functions.convert_to_units(temperature, KELVIN)
    def set_temperature(self, temperature):
        self.temperature = functions.convert_to_units(temperature, KELVIN)
    def set_pressure(self, pressure):
        self.pressure = functions.convert_to_units(pressure, KELVIN)
    def set_restraints(self, restraint):
        if restraint not in bss.Protocol.Equilibration.restraints():
            raise IOError(f"{restraint} is not in supported list of restraints {bss.Protocol.Equilibration.restraints()}")
        else:
            self.restraints = restraint
    def set_configuration(self, config):
        self.configuration = config
    def set_name(self, name):
        self.name = name
    def get_minimisation_steps(self):
        return self.minimisation_steps
    def get_minimisation_stepsize(self):
        return self.minimisation_stepsize
    def get_minimisation_tolerance(self):
        return self.minimisation_tolerance
    def get_short_nvt(self):
        return self.short_nvt_runtime
    def get_nvt(self):
        return self.nvt_runtime
    def get_npt(self):
        return self.npt_runtime
    # def get_equilibration_dir(self):
    #     return self.equilibration_dir
    def get_start_temperature(self):
        return self.start_temperature
    def get_end_temperature(self):
        return self.end_temperature
    def get_temperature(self):
        return self.temperature
    def get_pressure(self):
        return self.pressure
    def get_name(self):
        return self.name

def main():
    pass


if __name__ == "__main__":
    main()