"""
Prepare ligands and protein for AFE calculations
"""

def prepare_meze(Network): 
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
    ligands_datfile = Network.create_ligand_dat_file()
    network_forward, network_backward = Network.create_network_files()
    protocol_file = Network.create_dat_file()
    return Network

def main():
    pass


if __name__ == "__main__":
    main()