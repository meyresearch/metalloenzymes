from Network import Network
import os
import BioSimSpace as bss
import numpy as np
import functions
from definitions import PICOSECOND
import csv
import MDAnalysis as mda
import MDAnalysis.analysis.distances 
import pandas as pd


class Meze(Network):

    def __init__(self, protein_file, metal="ZN", cut_off=2.6, force_constant_0=100,
                 water_file=None, workdir=os.getcwd(), qmmm_inputs=None, equilibration_path=None, output=None,
                 ligand_path=os.getcwd()+"/inputs/ligands/", ligand_charge=0, ligand_ff="gaff2",
                 group_name=None, protein_path=os.getcwd()+"/inputs/protein/", water_model="tip3p", protein_ff="ff14SB", 
                 engine=None, sampling_time=10, box_edges=20, box_shape="cubic", min_steps=1000, short_nvt=0, nvt=1, npt=1, 
                 min_dt=None, min_tol=None, repeats=0, temperature=300, pressure=1, threshold=None, n_normal=None, n_difficult=None):
        
        super().__init__(workdir=workdir, ligand_path=ligand_path, group_name=group_name, protein_file=protein_file, protein_path=protein_path, 
                         water_model=water_model, ligand_ff=ligand_ff, protein_ff=protein_ff, ligand_charge=ligand_charge, equilibration_path=equilibration_path, outputs=output,
                         engine=engine, sampling_time=sampling_time, box_edges=box_edges, box_shape=box_shape, min_steps=min_steps, short_nvt=short_nvt, nvt=nvt, npt=npt, 
                         min_dt=min_dt, min_tol=min_tol, repeats=repeats, temperature=temperature, pressure=pressure, threshold=threshold, n_normal=n_normal, n_difficult=n_difficult)
        
        self.md_time = functions.convert_to_units(sampling_time, PICOSECOND)
        self.universe = mda.Universe(self.protein_file, format="pdb")
        self.force_constant_0 = functions.check_float(force_constant_0)

        if qmmm_inputs:
            self.qmmm_input_directory = functions.path_exists(qmmm_inputs)
        else:
            os.rmdir(self.afe_input_directory)
            self.qmmm_input_directory = self.create_directory(f"/qmmm_input_files/")
        
        if output:
            self.output_directory = functions.path_exists(output)
        else:
            self.output_directory = self.create_directory(f"/outputs/")
            
        if water_file:
            self.xtal_water = bss.IO.readMolecules(water_file)

        if metal.upper() == "ZN":
            self.metal_resname = metal
            self.metals = self.universe.select_atoms(f"resname {self.metal_resname}")
            self.metal_resids = self.metals.resids
            self.metal_atomids = self.metals.atoms.ids
        else:
            print(f"Metal {metal} is not supported yet.")
        
        self.cut_off = functions.check_positive(functions.check_float(cut_off))


    def prepare_meze(self):
        """
        Prepare QM/MM calculations by creating ligand and protocol dat files.

        Parameters:
        -----------

        Return:
        -------
        self: Network
            (prepared) Network object
        """
        self.protein_water_complex = self.protein.create_complex()
        self.prepared_protein = self.protein.tleap(self.protein_water_complex)
        self.log_directory = self.create_directory("/logs/")
        self.protocol_file = self.create_protocol_file()
        self.prepared = True
        return self


    def create_protocol_file(self):
            """
            Create protocol.dat file for QM/MM runs

            Parameters:
            -----------
            Network: Network
                Network class object

            Return:
            -------
            protocol_file: str
                protocol datafile
            """
            protocol = [f"metal = {self.metal_resname}",
                        f"cut-off = {self.cut_off}",
                        f"force constant 0 = {self.force_constant_0}",
                        f"group name = {self.group_name}",
                        f"ligand forcefield = {self.ligand_forcefield}", 
                        f"ligand charge = {self.ligand_charge}",
                        f"prepared protein file = {self.prepared_protein}",
                        f"protein input file = {self.protein_file}",
                        f"protein forcefield = {self.protein_forcefield}", 
                        f"water model = {self.water_model}", 
                        f"box edges = {self.box_edges}", # in angstrom 
                        f"box shape = {self.box_shape}", 
                        f"minimisation steps = {self.min_steps}",
                        f"minimisation stepsize = {self.min_dt}",
                        f"minimisation tolerance = {self.min_tol}",
                        f"short nvt = {self.short_nvt._value}",
                        f"nvt = {self.nvt._value}",
                        f"npt = {self.npt._value}",
                        f"temperature = {self.temperature._value}",
                        f"pressure = {self.pressure._value}",
                        f"sampling time = {self.md_time._value}",
                        f"engine = {self.md_engine}",
                        f"outputs = {self.output_directory}",
                        f"repeats = {self.n_repeats}",
                        f"project directory = {self.workding_directory}",
                        f"equilibration directory = {self.equilibration_directory}",
                        f"ligand directory = {self.ligand_path}",
                        f"protein directory = {self.protein_path}",
                        f"log directory = {self.log_directory}",
                        f"qmmm input directory = {self.qmmm_input_directory}"]

            protocol_file = self.qmmm_input_directory + "/protocol.dat"
            with open(protocol_file, "w") as file:
                writer = csv.writer(file)
                for protocol_line in protocol:
                    writer.writerow([protocol_line])
            return protocol_file
    

    def set_universe(self, file_name):
        """
        Take a file name and create an MDAnalysis Universe from prm7 and rst7 files.
        Set the universe attribute.

        Parameters:
        -----------
        file_name: str
            file name without extension

        Return:
        -------
        self: Meze.Meze()
            Meze class object
        """
        self.universe = mda.Universe(file_name + ".prm7", file_name + ".rst7",
                                     topology_format="PARM7", format="RESTRT")
        self.metals = self.universe.select_atoms(f"resname {self.metal_resname}")
        self.metal_resids = self.metals.resids
        self.metal_atomids = self.metals.atoms.ids
        return self
    

    def get_active_site(self, ligand_name):
        """
        Find active site residues defined as around metal ions defined by the given cutoff value.

        Parameters:
        -----------
        ligand_name: str
            name of ligand

        Return:
        -------
        mda.AtomGroup
            active site atoms
        """
        selection = f"resname {self.metal_resname} or around {self.cut_off} (resname {self.metal_resname})"
        return self.universe.select_atoms(selection)


    def get_qm_region(self, active_site):
        """
        Take an atom group of the active site and return a dictionary of the QM region selections.

        Parameters:
        -----------
        active_site: mda.AtomGroup
            group of atoms contained in the active site

        Return:
        -------
        qm_region: dict
            qm region selection for QM/MM
        """
        protein = self.universe.select_atoms("protein")
        qm_region = {}
        whole_residues = []
        atom_ids = []
        for residue in active_site.residues:
            if residue in protein.residues:
                selection = get_side_chain_selection(self.universe, residue)
                atom_ids.append(selection)
            else:
                whole_residues.append(str(residue.resid))
        qm_region["whole_residues"] = whole_residues
        qm_region["atom_ids"] = atom_ids
        return qm_region
    

    def dftb3(self, qm_region):
        """
        Create QM/MM options for DFTB3 to be used in QM/MM input files

        Parameters:
        -----------
        qm_region: dict
            dictionary of the qm region selections

        Return:
        -------
        options: dict
            DFTB3 options
        """
        whole_residues = ",".join(qm_region["whole_residues"])
        atom_ids = ",".join(qm_region["atom_ids"])
        options = {"qmmask": f"':{whole_residues}|(@{atom_ids})'\n",
                   "writepdb": "1\n",
                   "qmcharge": str(self.ligand_charge) + "\n",
                   "qm_theory": "'DFTB3'\n",
                   "qmshake": "0\n",
                   "qm_ewald": "1",
                   "qm_pme": "1"}
        return options
    
 
    def get_metal_ligands(self):
        """
        Create a dictionary of metal ion(s) and their coordinating ligand(s)

        Parameters:
        -----------

        Return:
        -------
        metal_ligands: dict
            dictionary where keys are metal ions (differentiated by subscripts) and values are the atom groups containing coordinating ligands for that ion
        """
        metal_ligands = {}
        for i in range(len(self.metal_resids)):
            ligands = self.universe.select_atoms(f"not (name {self.metal_resname} or element H) and sphzone {self.cut_off} (resid {self.metal_resids[i]})")
            key = self.metal_atomids[i] 
            metal_ligands[key] = ligands
        for key, value in metal_ligands.items():
            print(f"{key}: {value}")
        return metal_ligands


    def model_0(self, active_site):
        pass


    def write_restraints_file_0(self):
        """
        Write an Amber-compatible restraint file for model 0. 

        Parameters:
        -----------

        Return:
        -------
        """
        metal_ligands = self.get_metal_ligands()
        protein = self.universe.select_atoms("protein")
        atom_ids = []
        r1, r2, r3, r4 = [], [], [], []
        for i in range(len(self.metal_atomids)):
            metal_id = self.metal_atomids[i]
            for ligating_atom in metal_ligands[metal_id]:
                if ligating_atom in protein or ligating_atom.resname == "WAT":

                    atom_ids.append(f"iat={metal_id},{ligating_atom.id}")

                    atom_group_1 = self.universe.select_atoms(f"resid {self.metal_resids[i]}")
                    atom_group_2 = self.universe.select_atoms(f"resid {ligating_atom.resid} and name {ligating_atom.name}")
                    distance = MDAnalysis.analysis.distances.dist(atom_group_1, atom_group_2)[-1][0]
                    linear_response_region_lower_bound = np.round(distance - 1.0, decimals=2)
                    flat_region_lower_bound = np.round(distance - 0.5, decimals=2)
                    flat_region_upper_bound = np.round(distance + 0.5, decimals=2)
                    linear_response_region_upper_bound = np.round(distance + 1.0, decimals=2)
                    r1.append("r1="+str(linear_response_region_lower_bound))
                    r2.append("r1="+str(flat_region_lower_bound))
                    r3.append("r1="+str(flat_region_upper_bound))
                    r4.append("r1="+str(linear_response_region_upper_bound))

        rk2 = ["rk2="+str(self.force_constant_0)] * len(atom_ids)
        rk3 = ["rk3="+str(self.force_constant_0)] * len(atom_ids)

        with open(self.qmmm_input_directory + "restraints_0.RST", "w") as file:
            file.write(f"# Harmonic bond restraints between {self.metal_resname} and coordinating protein residues\n")
            for i in range(len(atom_ids)):
                line = f"&rst {atom_ids[i]}, {r1[i]}, {r2[i]}, {r3[i]}, {r4[i]}, {rk2[i]}, {rk3[i]},/\n"
                file.write(line)



    
def get_side_chain_selection(universe, residue):
    """
    Take a residue and return the atom indices of the side chain of that residue.

    Parameters:
    -----------
    universe: mda.Universe
    residue: mda.residue

    Return:
    -------
    str 
        atom selection in the format "{first_atom}-{last_atom}"
    """
    n_terminus = "name N or name H"
    alpha_carbon = "name CA or name HA"
    c_terminus = "name C or name O"
    atoms_in_residue = universe.select_atoms(f"resid {residue.resid}")
    qm_region_for_residue = list(atoms_in_residue.select_atoms(f"not ({n_terminus} or {alpha_carbon} or {c_terminus})").ids)
    first_atom = qm_region_for_residue[0]
    last_atom = qm_region_for_residue[-1]
    return f"{first_atom}-{last_atom}"
