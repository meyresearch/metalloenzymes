from Network import Network
import os
import BioSimSpace as bss
import Ligand
import functions
import multiprocessing
import tqdm


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
 