from Network import Network
import os
import BioSimSpace as bss
import Ligand
import functions


class Meze(Network):

    def __init__(self, workdir, ligand_path, group_name, protein_file, protein_path, water_model, ligand_ff, protein_ff, ligand_charge, 
                 sampling_time, box_edges, box_shape, min_steps, short_nvt, nvt, npt, 
                 min_dt, min_tol, water_file, repeats=0, temperature=300, pressure=1, threshold=None, n_normal=None, n_difficult=None, engine=None):
        
        super().__init__(workdir, ligand_path, group_name, protein_file, protein_path, water_model, ligand_ff, protein_ff, ligand_charge, 
                         engine, sampling_time, box_edges, box_shape, min_steps, short_nvt, nvt, npt, 
                         min_dt, min_tol, repeats, temperature, pressure, threshold, n_normal, n_difficult)
        
        os.rmdir(self.afe_input_directory)
        self.qmmm_input_directory = self.create_directory(f"/qmmm_inputs/")
        self.output_directory = self.create_directory(f"/outputs/")
        self.water_file = water_file
        self.xtal_water = bss.IO.readMolecules(self.water_file)
        self.bound_ligands = [self.create_bound_complex(i) for i in range(self.n_ligands)]
        

    def create_bound_complex(self, index):
        """
        Combine protein, ligand and crystal waters

        Parameters:
        -----------
        index: int
            Ligand index for sorting through Ligands
        
        Return:
        -------
        Ligand.Ligand
            bound ligand object with the crystal waters
        """        
        ligand_number = self.names[index].split("_")[-1]
        complex = self.protein.molecule + self.ligand_molecules[index] + self.xtal_water # ANTECHAMBER THE LIGS FIRST?
        file = bss.IO.saveMolecules(self.protein_path + "complex_" + ligand_number, complex, fileformat="pdb")[0]
        return Ligand.Ligand(file=file)

        
    def solvate_bound(self, index):
        """
        Solvate bound systems.

        Parameters:
        -----------
        index: int
            Ligand indices for sorting through Network.names and Network.ligands
        Return:
        -----
        solvated_system: Ligand
            (solvated) Ligand object whose file attribute is the prm7 and rst7 files         
        """

        ligand_number = self.names[index].split("_")[-1]
        print(f"Solvating bound ligand {ligand_number}")
        
        prepared_protein = self.protein.tleap(self.bound_ligands[index].file)

        # ligand_parameters = ligand.parameterise(self.ligand_forcefield, self.ligand_charge)     
        system_parameters = self.protein.get_prepared_protein(prepared_protein)
        bound_box, bound_box_angles = self.create_box(system_parameters)
        solvated_molecules = bss.Solvent.solvate(model=self.protein.water_model,
                                                 molecule=system_parameters,
                                                 box=bound_box,
                                                 angles=bound_box_angles)
        
        system_savename = self.protein_path + "system_" + ligand_number + "_solvated"
        solvated_files = bss.IO.saveMolecules(system_savename, solvated_molecules, ["PRM7", "RST7"])
        solvated_system = Ligand.Ligand(file=solvated_files, parameterised=True)
        return solvated_system