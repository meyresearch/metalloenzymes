"""
Ligands class object
"""
import functions
import BioSimSpace as bss
import Network 
import csv


class Ligand(object):
    """
    Ligand object for AFE preparation.

    Parameters:
    -----------
    path: str 
        path to /inputs/ligands
    forcefield: str
        ligand forcefield
 
    Methods:
    -------

    """
    def __init__(self, file, parameterised=False):
        """
        Class constructor
        """
        self.parameterised = parameterised
        if not self.parameterised:
            self.file = functions.file_exists(file)
        elif self.parameterised:
            self.file = file # file arg will be a list of prm7, rst7
        self.molecule = self.get_ligand()
        self.name = self.get_name()


    def get_ligand(self):
        """
        Get ligand molecule in path as bss.Molecule object.

        Return:
        -------
        ligands: bss.Molecule
            ligand molecule
        """
        ligand = bss.IO.readMolecules(self.file)[0]
        return ligand


    def get_name(self):
        """
        Get ligand names in path.

        Return:
        -------
        names: list
        """
        if self.parameterised:
            names = functions.get_filenames(self.file[0])
        else:
            names = functions.get_filenames(self.file)
        return names
    

    def parameterise(self, forcefield, charge):
        self.parameters = None
        if forcefield == "gaff":
            self.parameters = bss.Parameters.gaff(molecule=self.molecule, net_charge=charge).getMolecule()
        elif forcefield == "gaff2":
            self.parameters = bss.Parameters.gaff2(molecule=self.molecule, net_charge=charge).getMolecule()
        else:
            print(f"Forcefield {forcefield} is not supported for ligands.")
        return self.parameters  


    def get_system(self):
        """
        Get ligand system as a BSS system object

        Return:
        -------
        ligands: bss.System
            ligand system
        """
        system = bss.IO.readMolecules(self.file)
        return system
    

def main():
    lig = Ligand("/home/jguven/projects/alchemistry/add_caps_to_kpc2/inputs/ligands/ligand_1.sdf")
    # ntwrk = Network.create_network(mols, names, lig.path)
    print(lig.molecule, lig.name)

if __name__ == "__main__":
    main()