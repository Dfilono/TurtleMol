import math

def AtomFillSphere(sphere, og, radii, max_attempts, tol):
    filled = []
    attempts = 0

    while attempts < max_attempts:
        new_mol = []

        for i, atom in enumerate(og):
            # Calculate the relative position of each atom within the molecule
            rel_x, rel_y, rel_z = og[i][1:]
            atom_type = atom[0]

            # Calculate the coordiantes within the sphere
            x = sphere.x + rel_x + tol
            y = sphere.y + rel_y + tol
            z = sphere.z + rel_z + tol

            # Adjust for atomic radii
            atom_radius = radii.get(atom_type, 0.0) # Get the radius for the atom type

            # Check if the new atom fits within the sphere
            if sphere.contains_points(x, y, z, atom_radius):
                new_atom = (atom_type, x, y, z)
                print(new_atom)

            if not is_overlap_atom(new_atom, filled, radii, tol):
                new_mol.append(new_atom)
        
        if len(new_mol) == len(og):
            filled.append(new_mol)
            attempts += 1

    return filled

def AtomDefiniteSphere():
    pass

def MoleculeFillSphere(sphere, og, radii, max_attempts, tol):
    filled = []
    attempts = 0

    while attempts < max_attempts:
        new_mol = []

        for i, atom in enumerate(og):
            # Calculate the relative position of each atom within the molecule
            rel_x = atom[1]
            rel_y = atom[2]
            rel_z = atom[3]
            atom_type = atom[0]

            # Calculate the coordiantes within the sphere
            x = sphere.x + rel_x
            y = sphere.y + rel_y
            z = sphere.z + rel_z

            # Adjust for atomic radii
            atom_radius = radii.get(atom_type, 0.0) # Get the radius for the atom type

            # Check if the new atom fits within the sphere
            if (
                sphere.contains_points(x, y, z, atom_radius) and not
                is_overlap_molecule((atom_type, x, y, z), filled, radii, tol)
            ):
                new_atom = (atom_type, x, y, z)
                new_mol.append(new_atom)
            else:
                break # If any atom doesn't fit, discard the whol molecule
        
        if len(new_mol) == len(og):
            filled.append(new_mol)
            attempts += 1

    return filled

def MoleculeDefiniteSphere():
    pass

def is_overlap_molecule(new_atom, filled_atoms, radii, tol):
    print(new_atom)
    for molecule in filled_atoms:
        for atom1, atom2 in zip(new_atom, molecule):
            distance = math.sqrt(
                (atom1[1] - atom2[1])**2 +
                (atom1[2] - atom2[2])**2 +
                (atom1[3] - atom2[3])**2
            )

            if distance < (radii[atom1[0]] + radii[atom2[0]]) + tol:
                #print(f'Fail: dist = {distance}, Radius = {radii[new_atom[0]] + radii[atom[0]]}')
                return True
        #print(f'Pass: dist = {distance}, Radius = {radii[new_atom[0]] + radii[atom[0]]}')
    return False

def is_overlap_atom(new_atom, filled_atoms, radii, tol):
    for atom in filled_atoms:
        distance = math.sqrt(
            (new_atom[1] - atom[1])**2 +
            (new_atom[2] - atom[2])**2 +
            (new_atom[3] - atom[3])**2
        )

        if distance < (radii[new_atom[0]] + radii[atom[0]]) + tol:
            #print(f'Fail: dist = {distance}, Radius = {radii[new_atom[0]] + radii[atom[0]]}')
            return True
        #print(f'Pass: dist = {distance}, Radius = {radii[new_atom[0]] + radii[atom[0]]}')
    return False
