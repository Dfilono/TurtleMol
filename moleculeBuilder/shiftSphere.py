'''Module places molecules in a sphere'''

import random
from isOverlap import isOverlapAtom, isOverlapMolecule

def AtomFillSphere(numshifts, sphere, og, radii, tol):
    '''Fills sphere with single atoms'''
    filled = []

    for zshifts in range(numshifts):
        for yshifts in range(numshifts):
            for xshifts in range(numshifts):
                for atom in og:
                    # Calculate the relative position of each atom within the molecule
                    rel_x = atom[1]
                    rel_y = atom[2]
                    rel_z = atom[3]
                    atom_type = atom[0]

                    # Calculate the coordiantes within the sphere
                    x = (sphere.x - sphere.radius) + rel_x + xshifts
                    y = (sphere.y - sphere.radius) + rel_y + yshifts
                    z = (sphere.z - sphere.radius) + rel_z + zshifts

                    # Adjust for atomic radii
                    atom_radius = radii.get(atom_type, 0.0) # Get the radius for the atom type

                    # Check if the new atom fits within the sphere
                    if sphere.contains_points(x, y, z, atom_radius):
                        new_atom = (atom_type, x, y, z)

                        if not isOverlapAtom(new_atom, filled, radii, tol):
                            filled.append(new_atom)

    return filled

def AtomDefiniteSphere(numMol, maxattempts, og, sphere, radii, tol):
    '''Places a defined number of atoms in a sphere'''
    filled = []
    attempts = 0

    while len(filled) < numMol and attempts <= maxattempts:
        for atom in og:

            # Calculate the shift for each tile and new point
            new_x = atom[1] + random.uniform((sphere.x - sphere.radius), (sphere.x + sphere.radius))
            new_y = atom[2] + random.uniform((sphere.y - sphere.radius), (sphere.y + sphere.radius))
            new_z = atom[3] + random.uniform((sphere.z - sphere.radius), (sphere.z + sphere.radius))
            atom_type = atom[0]

            # Adjust for atomic radii
            atom_radius = radii.get(atom_type, 0.0) # Get the radius for the atom type

            # Check if the new atom fits within the sphere
            if sphere.contains_points(new_x, new_y, new_z, atom_radius):
                new_atom = (atom_type, new_x, new_y, new_z)

                if not isOverlapAtom(new_atom, filled, radii, tol):
                    filled.append(new_atom)

        attempts += 1

    return list(filled)

def MoleculeFillSphere(numshifts, sphere, og, radii, tol):
    '''Fills sphere with molecules'''
    filled = []

    for zshifts in range(numshifts):
        for yshifts in range(numshifts):
            for xshifts in range(numshifts):
                new_mol = []

                for i, atom in enumerate(og):
                    # Calculate the relative position of each atom within the molecule
                    rel_x = atom[1]
                    rel_y = atom[2]
                    rel_z = atom[3]
                    atom_type = atom[0]

                    # Calculate the coordiantes within the sphere
                    x = (sphere.x - sphere.radius) + rel_x + xshifts
                    y = (sphere.y - sphere.radius) + rel_y + yshifts
                    z = (sphere.z - sphere.radius) + rel_z + zshifts

                    # Adjust for atomic radii
                    atom_radius = radii.get(atom_type, 0.0) # Get the radius for the atom type

                    # Check if the new atom fits within the sphere
                    if sphere.contains_points(x, y, z, atom_radius):
                        new_atom = (atom_type, x, y, z)
                        new_mol.append(new_atom)
                    else:
                        break # If any atom doesn't fit, discard the whol molecule

                if not isOverlapMolecule(new_mol, filled, radii, tol):
                    if len(new_mol) == len(og):
                        filled.append(new_mol)

    return filled

def MoleculeDefiniteSphere(numMol, maxattempts, og, sphere, radii, tol):
    '''Places a defined number of molecules in a sphere'''
    filled = []
    attempts = 0

    while len(filled) < numMol and attempts <= maxattempts:
        new_mol = []
        shiftx = random.uniform((sphere.x - sphere.radius), (sphere.x + sphere.radius))
        shifty = random.uniform((sphere.y - sphere.radius), (sphere.y + sphere.radius))
        shiftz = random.uniform((sphere.z - sphere.radius), (sphere.z + sphere.radius))

        for atom in og:

            # Calculate the shift for each tile and new point
            new_x = atom[1] + shiftx
            new_y = atom[2] + shifty
            new_z = atom[3] + shiftz
            atom_type = atom[0]

            # Adjust for atomic radii
            atom_radius = radii.get(atom_type, 0.0) # Get the radius for the atom type
            # Check if the new atom fits within the sphere
            if sphere.contains_points(new_x, new_y, new_z, atom_radius):
                new_atom = (atom_type, new_x, new_y, new_z)
                new_mol.append(new_atom)
            else:
                break # If any atom doesn't fit, discard the whol molecule

        if not isOverlapMolecule(new_mol, filled, radii, tol):
            if len(new_mol) == len(og):
                filled.append(new_mol)

        attempts += 1

    return list(filled)
