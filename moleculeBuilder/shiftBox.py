import random
import math

def AtomsFillBox(x, y, z, tol, og, box, radii):
    filledAtom = []
    filledPositions = set()
    for zshift in range(z):
        for yshift in range(y):
            for xshift in range(x):
                for atom in og:
                    # Calculate the shift for each tile and new point
                    new_x = box.xCoord + atom[1] + xshift
                    new_y = box.yCoord + atom[2] + yshift
                    new_z = box.zCoord + atom[3] + zshift

                    # Adjust for atomic radiss
                    new_x_min = new_x - radii[atom[0]]
                    new_y_min = new_y - radii[atom[0]]
                    new_z_min = new_z - radii[atom[0]]
                    new_x_max = new_x + radii[atom[0]]
                    new_y_max = new_y + radii[atom[0]]
                    new_z_max = new_z + radii[atom[0]]

                    # Check if the new atom fits within the box
                    if  inBox(new_x_min, new_x_max, new_y_min, new_y_max, new_z_min, new_z_max, box):
                        new_atom = (atom[0], new_x, new_y, new_z)

                        if not is_overlap_atom(new_atom, filledAtom, radii, tol):
                            filledAtom.append(new_atom)
                            filledPositions.add(new_atom)

    return filledAtom

def AtomsDefiniteBox(numMol, maxattempts, og, box, radii, tol): # NOTE only works by randomly assigning atoms
    filledAtom = []
    filledPositions = set()

    attempts = 0

    while len(filledAtom) < numMol and attempts < maxattempts:
        attempts += 1
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
            if  inBox(new_x_min, new_x_max, new_y_min, new_y_max, new_z_min, new_z_max, box):
                new_atom = (atom[0], new_x, new_y, new_z)

                if not is_overlap_atom(new_atom, filledAtom, radii, tol):
                    filledAtom.append(new_atom)
                    filledPositions.add(new_atom)


    return list(filledAtom)

def MoleculesFillBox(x, y, z, tol, og, box, radii): 
    filledAtom = []
    filledPositions = set()

    for zshift in range(z):
        for yshift in range(y):
            for xshift in range(x):
                new_mol = []
                any_atom_outside = False

                for i, atom in enumerate(og):
                    #print(atom)
                    # Calculate the shift for each tile and new point
                    new_x = box.xCoord + atom[1] + xshift
                    new_y = box.yCoord + atom[2] + yshift
                    new_z = box.zCoord + atom[3] + zshift

                    # Adjust for atomic radiss
                    new_x_min = new_x - radii[atom[0]]
                    new_y_min = new_y - radii[atom[0]]
                    new_z_min = new_z - radii[atom[0]]
                    new_x_max = new_x + radii[atom[0]]
                    new_y_max = new_y + radii[atom[0]]
                    new_z_max = new_z + radii[atom[0]]

                    # Check if the new atom fits within the box
                    if  inBox(new_x_min, new_x_max, new_y_min, new_y_max, new_z_min, new_z_max, box):
                        new_atom = (atom[0], new_x, new_y, new_z)
                        new_mol.append(new_atom)
                    else:
                        any_atom_outside = True

                if not is_overlap_molecule(new_mol, filledAtom, radii, tol) and not any_atom_outside:
                    filledAtom.append(new_mol)
                    print("Added molecule:", new_mol)

    return filledAtom

def MoleculesDefiniteBox(numMol, maxattempts, og, box, radii, tol): 
    filledAtom = []
    filledPositions = set()

    attempts = 0

    while len(filledAtom) < numMol and attempts < maxattempts:
        attempts += 1
        new_mol = []
        any_atom_outside = False

        xshift = random.uniform(0, box.length)
        yshift = random.uniform(0, box.width)
        zshift = random.uniform(0, box.height)

        for atom in og:

            # Calculate the shift for each tile and new point
            new_x = atom[1] + xshift
            new_y = atom[2] + yshift
            new_z = atom[3] + zshift

            # Adjust for atomic radiss
            new_x_min = new_x - radii[atom[0]]
            new_y_min = new_y - radii[atom[0]]
            new_z_min = new_z - radii[atom[0]]
            new_x_max = new_x + radii[atom[0]]
            new_y_max = new_y + radii[atom[0]]
            new_z_max = new_z + radii[atom[0]]

            # Check if the new atom fits within the box
            if  inBox(new_x_min, new_x_max, new_y_min, new_y_max, new_z_min, new_z_max, box):
                new_atom = (atom[0], new_x, new_y, new_z)
                new_mol.append(new_atom)
            else:
                any_atom_outside = True

        if not is_overlap_molecule(new_mol, filledAtom, radii, tol) and not any_atom_outside:
            filledAtom.append(new_mol)
            print("Added molecule:", new_mol)

    return list(filledAtom)

def is_overlap_molecule(new_atom, filled_atoms, radii, tol):
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

def inBox(xmin, xmax, ymin, ymax, zmin, zmax, box):
    if (
        xmin >= box.xCoord and xmax <= box.xCoord + box.length and
        ymin >= box.yCoord and ymax <= box.yCoord + box.width and
        zmin >= box.zCoord and zmax <= box.zCoord + box.height
    ):
        return True
    else:
        return False