"""
Alchemical Free Energy class for free energy simulations
"""
import pathlib
import csv

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
    def __init__(self, path, engine, sampling_time, box_edges, box_shape):
        self.path = path
        self.engine = engine
        self.time = sampling_time
        self.box = (box_edges, box_shape)

    def create_directory(self):
        """
        Create AFE working directory in path.
        Set home attribute.

        Parameters:
        -----------
        Return:
        -------
        directory: str
            full path to afe directory
        """
        try:
            directory = self.path + "/afe/"
            pathlib.Path(directory).mkdir(parents=False, exist_ok=False)
        except FileNotFoundError as e:
            print(f"Could not create afe directory. Pathlib raised error: {e}")
        except FileExistsError as e:
            print(f"Could not create afe directory. Pathlib raised error: {e}")
        self.home = directory
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
                    f"box edges = {self.box[0]}*angstrom", 
                    f"box shape = {self.box[1]}", 
                    f"protocol = default",
                    f"sampling = {self.time}*ns",
                    f"engine = {self.engine}"]
        self.protocol_file = self.home + "/protocol.dat"

        with open(self.protocol_file, "w") as file:
            writer = csv.writer(file)
            for protocol_line in protocol:
                writer.writerow([protocol_line])
 
    

def main():
    pass


if __name__ == "__main__":
    main()