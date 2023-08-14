"""
Prepare ligands and protein for AFE calculations
"""

def prepare_meze(Protein, Network, AFE): 
    """
    Prepare ligands and protein for AFE calculations.

    Parameters:
    -----------
    Protein: Protein
        Protein class object
    Ligands: Ligands
        Ligands class object
    AFE: AlchemicalFreeEnergy
        AlchemicalFreeEnergy class object

    Return:
    -------
    """
    network_dictionary = Network.create_dictionary()
    #TODO Edit dictionary?
    protein_water_complex_file = Protein.create_complex()
    prepared_protein_file = Protein.tleap(protein_water_complex_file)

    ligands_datfile = Network.create_ligand_dat_file(AFE.afe_dir)
    network_forward, network_backward = Network.create_network_files(AFE.engine)
    protocol_file = AFE.create_dat_file(Protein=Protein, Network=Network)


def main():
    pass


if __name__ == "__main__":
    main()