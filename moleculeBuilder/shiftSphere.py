import math
import random

def AtomFillSphere(numshifts, sphere, og, radii, tol):
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

                        if not is_overlap_atom(new_atom, filled, radii, tol):
                            filled.append(new_atom)

    return filled

def AtomDefiniteSphere(numMol, maxattempts, og, sphere, radii, tol):
    filled = []
    attempts = 0

    while len(filled) < numMol and attempts <= maxattempts:
        for atom in og:

            # Calculate the shift for each tile and new point
            new_x = atom[1] + random.uniform((sphere.x - sphere.radius), sphere.radius)
            new_y = atom[2] + random.uniform((sphere.y - sphere.radius), sphere.radius)
            new_z = atom[3] + random.uniform((sphere.z - sphere.radius), sphere.radius)
            atom_type = atom[0]

            

            # Adjust for atomic radii
            atom_radius = radii.get(atom_type, 0.0) # Get the radius for the atom typ
            # Check if the new atom fits within the sphere
            if sphere.contains_points(new_x, new_y, new_z, atom_radius):
                new_atom = (atom_type, new_x, new_y, new_z)

                if not is_overlap_atom(new_atom, filled, radii, tol):
                    filled.append(new_atom)
            
        attempts += 1

    return list(filled)

def MoleculeFillSphere(numshifts, sphere, og, radii, tol):
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

                if not is_overlap_molecule(new_mol, filled, radii, tol):
                    if len(new_mol) == len(og):
                        filled.append(new_mol)

    return filled

def MoleculeDefiniteSphere(numMol, maxattempts, og, sphere, radii, tol):
    filled = []
    attempts = 0

    while len(filled) < numMol and attempts <= maxattempts:
        new_mol = []
        shiftx = random.uniform((sphere.x - sphere.radius), sphere.radius)
        shifty = random.uniform((sphere.y - sphere.radius), sphere.radius)
        shiftz = random.uniform((sphere.z - sphere.radius), sphere.radius)

        for atom in og:

            # Calculate the shift for each tile and new point
            new_x = atom[1] + shiftx
            new_y = atom[2] + shifty
            new_z = atom[3] + shiftz
            atom_type = atom[0]            

            # Adjust for atomic radii
            atom_radius = radii.get(atom_type, 0.0) # Get the radius for the atom typ
            # Check if the new atom fits within the sphere
            if sphere.contains_points(new_x, new_y, new_z, atom_radius):
                new_atom = (atom_type, new_x, new_y, new_z)
                new_mol.append(new_atom)
            else:
                break # If any atom doesn't fit, discard the whol molecule
        
        if not is_overlap_molecule(new_mol, filled, radii, tol):
            if len(new_mol) == len(og):
                filled.append(new_mol)

        attempts += 1

    return list(filled)

def is_overlap_molecule(new_atom, filled_atoms, radii, tol):
    #print(new_atom)
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
