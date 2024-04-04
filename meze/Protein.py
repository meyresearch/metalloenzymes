"""
Protein class object
"""
import os
import functions
import BioSimSpace as bss
import argparse


def get_water_file(working_directory, name="water.pdb"):
    """
    Get water.pdb file from current working directory and check it exists
    #TODO maybe get it as input at least the name of the file

    Parameters:
    -----------
    Return:
    -------
    water_file: str or None
        crystal waters in pdb file or None if water file doesn't exist in working directory
    """
    water_file = working_directory + "/" + name
    try:
        water_file = functions.file_exists(water_file)
    except argparse.ArgumentTypeError:
        print("\n")
        print("No separate water file detected in working directory. Continuing without it.")
        water_file = None
    return water_file


class Protein(object):
    """
    Protein object for AFE preparation

    Attributes:
    -----------
    name: str
        name of the protein: used for saving meze output files
    protein_file: str
        prepared protein file
    path: str
        path to protein prep directory
    forcefield: str
        protein forcefield, default is ff14SB
    water_model: str
        water model, default is tip3p
        
    Methods:
    -------
    create_complex() -> str:
        Combine protein file with crystallographic waters.
        Set complex attribute.
    write_tleap_input() -> str:
        Write tleap input file for setting up protein with given parameters
    tleap() -> None:
        Run tleap for water-protein complex.
    """
    def __init__(self, name, protein_file, path, forcefield, water_model, parameterised=False):
        """
        Class constructor
        """
        self.name = name
        self.path = functions.path_exists(path)
        self.parameterised = parameterised
        if not self.parameterised:
            self.file = functions.file_exists(protein_file)
        elif self.parameterised:
            self.file = protein_file 
        self.molecule = self.get_molecule()
        self.forcefield = forcefield
        self.water_model = water_model
        self.prepared = self.path + "/" + self.name + "_tleap"
    

    def get_molecule(self):
        """
        Get protein molecule(s)

        Parameters:
        -----------

        Return:
        -------
        bss.Molecules
            molecules from the protein file(s)
        """
        return bss.IO.readMolecules(self.file)

    
    def create_complex(self):
        """
        Combine protein file with crystallographic waters.
        Set complex attribute.

        Parameters:
        -----------
        Return:
        -------
        pdb_file: str
            pdb file with water and protein combined
        """
        water_file = get_water_file(self.path)
        output = self.path + "/" + self.name
        if water_file is not None:
            self.xtal_water = bss.IO.readMolecules(water_file)
            complex = self.molecule + self.xtal_water
            file = bss.IO.saveMolecules(output, complex, fileformat="pdb")[0]
        else:
            file = bss.IO.saveMolecules(output, self.molecule, fileformat="pdb")[0]
        
        return file


    def tleap(self, complex_file):
        """
        Run tleap for protein-water complex
        
        Parameters:
        -----------
        complex_file: str
            protein-water complex pdb file
        Return:
        -------
        save_file: str
            name of prepared protein-water system
        """
        tleap_input_file = self.path + "/tleap.in"
        tleap_output_file = self.path + "/tleap.out"
        #TODO Check if this works:
        # save_file = self.path + "/" + self.name
        save_file = functions.get_filename(complex_file)
        if not os.path.isfile(f"{save_file}_tleap.prm7") or not os.path.isfile(f"{save_file}_tleap.rst7"):
            with open(tleap_input_file, "w") as tleap_in:
                tleap_in.write(f"source leaprc.protein.{self.forcefield}\n")
                tleap_in.write(f"source leaprc.water.{self.water_model}\n")
                tleap_in.write(f"complex = loadpdb {complex_file}\n")
                tleap_in.write(f"saveamberparm complex {save_file}_tleap.prm7 {save_file}_tleap.rst7\n")
                tleap_in.write("quit")
            work_dir = os.getcwd()
            os.chdir(self.path)
            os.system(f"tleap -s -f {tleap_input_file} > {tleap_output_file}")
            os.chdir(work_dir)
        return self.path + save_file + "_tleap"
        

def main():
    prot = Protein(name="test", 
                   protein_file="/home/jguven/projects/alchemistry/add_caps_to_kpc2/inputs/protein/kpc2.input.pdb",
                   path="/home/jguven/projects/alchemistry/add_caps_to_kpc2/inputs/protein/",
                   forcefield="ff14SB",
                   water_model="tip3p")
    complex = prot.create_complex()
    command = prot.write_tleap_input(complex_file=complex)
    prot.tleap(command)


if __name__ == "__main__":
    main()
    
