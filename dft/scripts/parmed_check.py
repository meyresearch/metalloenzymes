def check_parmed_output(parmed_out: str) -> None:
    """
    Open a parmed output file and perform checks as described in http://ambermd.org/tutorials/advanced/tutorial20/mcpbpy.htm
    """
    with open(parmed_out, "r") as file:
        lines = file.readlines()
    for i in range(len(lines)):
        if "Reading" in lines[i]:
            start_line = i + 1

    info = [line.split() for line in lines[start_line:]]
    sections = []
    for i in range(len(info)):
        if len(info[i]) > 0:
            first_element = info[i][0]
            if "atom" in first_element.lower():
                sections.append(i)

    bonds = info[sections[0] + 1:sections[1]]
    force_constants, equilibrium_distances = [], []
    atom_info = []
    for bond_info in bonds:
        if len(bond_info) != 0:
            atoms = bond_info[0] + "-" + bond_info[1] + " " + bond_info[4] + "-" + bond_info[5]
            atom_info.append(atoms)
            equilibrium_distances.append(bond_info[-2])
            force_constants.append(bond_info[-1])

    for i in range(len(force_constants)):
        if float(equilibrium_distances[i]) >= 2.8:
            print(f"WARNING: The R eq. of {equilibrium_distances[i]} for bond {atom_info[i]} is greater than the expected 2.8 Å")
        if float(force_constants[i]) >= 200:
            print(f"WARNING: The force constant of {force_constants[i]} for bond {atom_info[i]} is greater than the expected 200 kcal / mol Å^2")

    angles = info[sections[1] + 1:sections[2]]
    angle_force_constants, equilibrium_angles = [], []
    atom_angle_info = []
    for angle_info in angles:
        if len(angle_info) != 0:
            atoms = angle_info[0] + "-" + angle_info[1] + " " + angle_info[4] + "-" + angle_info[5] + " " + angle_info[8] + "-" + angle_info[9]
            atom_angle_info.append(atoms)
            equilibrium_angles.append(angle_info[-1])
            angle_force_constants.append(angle_info[-2])

    for i in range(len(angle_force_constants)):
        if float(equilibrium_angles[i]) < 100:
            print(f"WARNING: The theta eq. of {equilibrium_angles[i]} for angle {atom_angle_info[i]} is less than the expected 100 rad.")
        if float(angle_force_constants[i]) > 100:
            print(f"WARNING: The force constant of {angle_force_constants[i]} for angle {atom_angle_info[i]} is greater than the expected 100 kcal / mol rad^2")

    dihedrals = info[sections[2] + 1:sections[3]]
    dihedral_atoms = []
    potential_barriers = []
    for dihedral in dihedrals:
        if len(dihedral) != 0 and "mask" not in dihedral:
            atoms = dihedral[0] + "-" + dihedral[1] + " " + dihedral[4] + "-" + dihedral[5] + " " + dihedral[8] + "-" + dihedral[9] + " " + dihedral[12] + "-" + dihedral[13]
            dihedral_atoms.append(atoms)
            potential_barriers.append(dihedral[16])

    for i in range(len(potential_barriers)):
        if float(potential_barriers[i]) > 0:
            print(f"WARNING: The height of the potential barrier of {potential_barriers[i]} for dihedral {dihedral_atoms[i]} is greater than the expected 0")

    resp_fitting = info[sections[3]+1:]
    resnames = []
    charges = []
    radii = []
    for fit in resp_fitting:
        if len(fit) > 1:
            resnames.append(fit[2])
            charges.append(fit[9])
            radii.append(fit[-2])

    for i in range(len(resnames)):
        if float(charges[i]) > 1:
            print(f"WARNING: The RESP charge of {charges[i]} for {resnames[i]} is greater than the expected value of +1")
        if float(radii[i]) < 1:
            print(f"WARNING: The LJ radius of {radii[i]} for {resnames[i]} is less than the expected value of 1 Å")


# ligands = [1, 2, 3, 4, 5, 6, 7, 9, 10, 11, 12, 13, 14, 15, 16]
# ligands = [1, 3, 4, 6, 7, 9, 10, 12, 13, 16]
ligands = [2, 5, 11, 14, 15]

for lig in ligands:
    print(f"Ligand {lig}, ZN1")
    parmed_file = f"../ligand_{lig}/parmed_1.out"
    check_parmed_output(parmed_file)
    print(f"Ligand {lig}, ZN2")
    parmed_file = f"../ligand_{lig}/parmed_2.out"
    check_parmed_output(parmed_file)
