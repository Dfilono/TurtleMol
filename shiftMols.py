import math
import random

def shiftMoleculesFill(x, y, z, tol, og, box, radii): # NOTE currently does not maintain relative position
    filledAtom = []
    filledPositions = set()

    relative_pos = [(atom[1] - og[0][1], atom[2] - og[0][2], atom[3] - og[0][3]) for atom in og]

    for zshift in range(z):
        for yshift in range(y):
            for xshift in range(x):
                new_mol = []
                any_atom_outside = False

                for i, atom in enumerate(og):
                    # Calculate the shift for each tile and new point
                    new_x = box.x + atom[1] + xshift * tol
                    new_y = box.y + atom[2] + yshift * tol
                    new_z = box.z + atom[3] + zshift * tol

                    # Adjust for atomic radiss
                    new_x_min = new_x - radii[atom[0]]
                    new_y_min = new_y - radii[atom[0]]
                    new_z_min = new_z - radii[atom[0]]
                    new_x_max = new_x + radii[atom[0]]
                    new_y_max = new_y + radii[atom[0]]
                    new_z_max = new_z + radii[atom[0]]

                    # Check if the new atom fits within the box
                    if  (
                        new_x_min >= box.x and new_x_max <= box.x + box.length and
                        new_y_min >= box.y and new_y_max <= box.y + box.width and
                        new_z_min >= box.z and new_z_max <= box.z + box.height
                    ):
                        new_atom = (atom[0], new_x, new_y, new_z)
                        new_mol.append(new_atom)
                    else:
                        any_atom_outside = True

                if not is_overlap(new_mol, filledAtom, radii, tol) and not any_atom_outside:
                    filledAtom.append(new_mol)
                    print("Added molecule:", new_mol)
                    #filledPositions.add(new_mol)

    return filledAtom

def shiftMoleculesDefinite(numMol, maxattempts, og, box, radii): # NOTE DOES NOT WORK
    filledAtom = []
    filledPositions = set()

    attempts = 0

    while len(filledAtom) < numMol and attempts < maxattempts:
        attempts += 1
        new_mol = []
        for atom in og:

            # Calculate the shift for each tile and new point
            new_x = atom[1] + random.uniform(0, box.length)
            new_y = atom[2] + random.uniform(0, box.width)
            new_z = atom[3] + random.uniform(0, box.height)

            # Adjust for atomic radiss
            new_x_min = new_x - radii[atom[0]]
            new_y_min = new_y - radii[atom[0]]
            new_z_min = new_z - radii[atom[0]]
            new_x_max = new_x + radii[atom[0]]
            new_y_max = new_y + radii[atom[0]]
            new_z_max = new_z + radii[atom[0]]

            # Check if the new atom fits within the box
            if  (
                new_x_min >= box.x and new_x_max <= box.x + box.length and
                new_y_min >= box.y and new_y_max <= box.y + box.width and
                new_z_min >= box.z and new_z_max <= box.z + box.height
            ):
                new_atom = (atom[0], new_x, new_y, new_z)
                new_mol.append(new_atom)

        if not is_overlap(new_mol, filledAtom, radii):
            filledAtom.append(new_mol)
            #filledPositions.add(new_mol)


    return list(filledAtom)

def is_overlap(new_atom, filled_atoms, radii, tol):
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