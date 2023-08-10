"""
Ligands class object
"""
import functions
import BioSimSpace as bss
import network

class Ligands(object):
    """
    Ligands object for AFE preparation.

    Parameters:
    -----------
    path: str 
        path to /inputs/ligands
    forcefield: str
        ligand forcefield
 
    Methods:
    -------

    """
    def __init__(self, path, forcefield):
        """
        Class constructor
        """
        self.path = functions.path_exists(path)
        self.forcefiedl = forcefield


    def get_files(self):
        """
        Read in all .sdf and .mol2 files in the provided path.
        Set files attribute.

        Parameters:
        -----------
        ligand_path: str
            path to directory containing mol2 and/or sdf files

        Return:
        -------
        list
            list of ligand filenames
        """
        # Adapted from dbmol.py:
        ligand_files = functions.read_files(f"{self.path}/*.sdf")
        ligand_files += functions.read_files(f"{self.path}/*.mol2")        
        if len(ligand_files) < 2:
            raise IOError(f"Path {self.path} must contain at least two sdf or mol2 files.")
        self.files = ligand_files
        return ligand_files

    
    def get_molecules(self):
        """
        Get ligand molecules in path as bss.Molecule objects.
        Set molecules attribute.

        Return:
        -------
        ligands: list
            ligand molecules
        """
        ligands = [bss.IO.readMolecules(file)[0] for file in self.files]
        self.molecules = ligands
        return ligands


    def get_names(self):
        """
        Get ligand names in path.
        Set names attribute.

        Return:
        -------
        names: list
        """
        names = [functions.get_filenames(file) for file in self.files]
        self.names = names
        return names


def main():
    ligs = Ligands(path="/home/jguven/projects/alchemistry/add_caps_to_kpc2/inputs/ligands",
                   forcefield="gaff2")
    files = ligs.get_files()
    mols = ligs.get_molecules()
    names = ligs.get_names()
    ntwrk = network.create_network(mols, names, ligs.path)


if __name__ == "__main__":
    main()