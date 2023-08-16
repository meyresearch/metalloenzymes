"""
Alchemical Free Energy class for free energy simulations
"""
import pathlib
import csv
import BioSimSpace as bss
from definitions import ANGSTROM
import functions 


class AlchemicalFreeEnergy(object):
    """
    _summary_

    Attributes:
    -----------
    object: 
        _description_

    Methods:
    -------
    """
    def __init__(self, path, engine, sampling_time, box_edges, box_shape, min, short_nvt, nvt, npt):
        self.path = functions.path_exists(path)
        self.engine = engine
        self.time = sampling_time
        self.box_shape = box_shape
        self.box_edges = box_edges
        self.min_steps = min
        self.short_nvt = functions.convert_to_units(short_nvt)
        self.nvt = functions.convert_to_units(nvt)
        self.npt = functions.convert_to_units(npt)
        self.afe_dir = self.create_directory("/afe/")
        self.equilibration_dir = self.create_directory("/equilibration/")
        

    def create_directory(self, name):
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
            directory = self.path + str(name)
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


def main():
    pass


if __name__ == "__main__":
    main()