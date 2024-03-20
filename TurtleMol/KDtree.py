import scipy.spatial
import numpy as np

def KDcheck(molecules, radii, tol=0.0):
    """
    Checks if there is any collision between atoms of different molecules based on atomic radii and a tolerance value.

    :param molecules: A list of molecules, each molecule is a list of atoms,
                      and each atom can be represented as [Name, X, Y, Z] or [Name, X, Y, Z, Other].
    :param radii: A dictionary where keys are atom names and values are their radii.
    :param tol: A user-defined minimum distance between any two molecules.
    :return: True if there is a collision, False otherwise.
    """
    allAtoms = []

    for mol_id, molecule in enumerate(molecules):
        for atom in molecule:
            if len(atom) >= 4:
                name, x, y, z = atom[:4]
                radius = radii.get(name, 0)  # Default to 0 if not found
                allAtoms.append((mol_id, name, np.array([x, y, z]), radius))
            else:
                raise ValueError("Atom data format not recognized")

    # Extract just the coordinates for building the KD-tree
    atomCoords = [atom[2] for atom in allAtoms]

    # Build the KD Tree
    tree = scipy.spatial.cKDTree(atomCoords)

    # Collision detection
    for i, (mol_id, atom_name, position, radius_i) in enumerate(allAtoms):
        # Find atoms within the sum of radii plus tolerance
        potential_collisions = tree.query_ball_point(position, r=radius_i + tol)

        for j in potential_collisions:
            if i == j:  # Skip the same atom
                continue

            other_mol_id, other_atom_name, _, radius_j = allAtoms[j]
            if mol_id == other_mol_id:  # Skip atoms from the same molecule
                continue

            # Calculate the effective collision distance including tolerance
            collision_distance = radius_i + radius_j + tol
            actual_distance = np.linalg.norm(position - allAtoms[j][2])

            if actual_distance < collision_distance:
                return True  # Collision detected

    return False  # No collision detected
