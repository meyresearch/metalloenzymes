"""
Protein class object
"""
import os
import functions
import BioSimSpace as bss


def get_water_file(working_directory, name="water.pdb"):
    """
    Get water.pdb file from current working directory and check it exists
    #TODO maybe get it as input at least the name of the file

    Parameters:
    -----------
    Return:
    -------
    water_pdb: str
        crystal waters in pdb file
    """
    water_file = working_directory + "/" + name
    return functions.file_exists(water_file)


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
    def __init__(self, name, protein_file, path, forcefield, water_model):
        """
        Class constructor
        """
        self.name = name
        self.path = functions.path_exists(path)
        self.molecule = bss.IO.readMolecules(functions.file_exists(protein_file))
        self.forcefield = forcefield
        self.water_model = water_model
        self.xtal_water = bss.IO.readMolecules(get_water_file(self.path))
    
    
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

        output = self.path + "/protein_water_complex"
        complex = self.molecule + self.xtal_water
        self.complex = complex
        bss.IO.saveMolecules(output, complex, fileformat="pdb")
        return output + ".pdb"


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
        save_file = self.path + "/" + self.name
        with open(tleap_input_file, "w") as tleap_in:
            tleap_in.write(f"source leaprc.protein.{self.forcefield}\n")
            tleap_in.write(f"source leaprc.water.{self.water_model}\n")
            tleap_in.write(f"complex = loadpdb {complex_file}\n")
            tleap_in.write(f"saveamberparm complex {save_file}_tleap.prm7 {save_file}_tleap.rst7\n")
            tleap_in.write("quit")
        os.system(f"tleap -s -f {tleap_input_file} > {tleap_output_file}")
        self.prepared = save_file + "_tleap"
        return self.prepared


    def get_prepared_protein(self):
        """
        Return a bss.Molecule object of the prepared protein-water complex

        Return:
        -------
        bss.Molecule: 
            prepared molecule 
        """
        return bss.IO.readMolecules([f"{self.prepared}.rst7", f"{self.prepared}.prm7"])[0]
        

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
    
