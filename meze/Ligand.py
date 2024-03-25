"""
Ligands class object
"""
import functions
import BioSimSpace as bss
import os


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
            names = functions.get_filename(self.file[0])
        else:
            names = functions.get_filename(self.file)
        return names
    

    def parameterise(self, forcefield, charge):
        """
        Parameterise ligand using BioSimSpace

        Parameters:
        -----------
        forcefield: 
            name of the ligand force field, current options are GAFF and GAFF2
        charge: 
            charge of the ligand

        Return:
        -------
        self.parameters: bss.Molecule
            parameterised ligand molecule object
        """
        self.parameters = None
        if forcefield == "gaff":
            self.parameters = bss.Parameters.gaff(molecule=self.molecule, net_charge=charge).getMolecule()
        elif forcefield == "gaff2":
            self.parameters = bss.Parameters.gaff2(molecule=self.molecule, net_charge=charge).getMolecule()
        else:
            print(f"Forcefield {forcefield} is not supported for ligands.")
        return self.parameters  


    def antechamber(self, charge, path, input_file=None, atom_type="gaff2", charge_method="bcc"):
        """
        Parameterise ligand using antechamber

        Parameters:
        -----------
        charge: float
            charge of the ligand

        Return:
        -------
        self: Ligand
            Ligand object
        """
        if not input_file:
            input_file = self.file
            
        else:
            input_file = functions.file_exists(input_file)

        input_file_extension = functions.get_file_extension(input_file)
        command = f"antechamber -fi {input_file_extension} -fo mol2 -i {input_file} -o {self.name}.mol2 -c {charge_method} -at {atom_type} -pf y -nc {charge}"
        work_dir = os.getcwd()
        os.chdir(path)
        os.system(command)
        os.chdir(work_dir)
        self.file = f"{path}/{self.name}.mol2"
        self.molecule = self.get_ligand()
        return self


    def parmcheck(self, path, atom_type="gaff2"):
        """
        Run parmchk2 from amber to create an frcmod file

        Parameters:
        -----------

        Return:
        -------
        str:
            full path to frcmod file
        """
        command = f"parmchk2 -i {self.name}.mol2 -o {self.name}.frcmod -f mol2 -s {atom_type}"
        work_dir = os.getcwd()
        os.chdir(path)
        os.system(command)
        os.chdir(work_dir)    
        return functions.get_files(f"{path}/{self.name}.frcmod")[0]


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
    print(lig.molecule, lig.name)

if __name__ == "__main__":
    main()