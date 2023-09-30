'''Module places molecules in a box'''

import random
from isOverlap import isOverlapAtom, isOverlapMolecule

def AtomsFillBox(x, y, z, tol, og, box, radii):
    '''Fills box with single atoms'''
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
                    if inBox(new_x_min, new_x_max, new_y_min, new_y_max, new_z_min, new_z_max, box):
                        new_atom = (atom[0], new_x, new_y, new_z)

                        if not isOverlapAtom(new_atom, filledAtom, radii, tol):
                            filledAtom.append(new_atom)
                            filledPositions.add(new_atom)

    return filledAtom

def AtomsDefiniteBox(numMol, maxattempts, og, box, radii, tol):
    '''Places a defined number of single atoms in box'''
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

                if not isOverlapAtom(new_atom, filledAtom, radii, tol):
                    filledAtom.append(new_atom)
                    filledPositions.add(new_atom)

    return list(filledAtom)

def MoleculesFillBox(x, y, z, tol, og, box, radii):
    '''Fills box with molecules'''
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

                if not isOverlapMolecule(new_mol, filledAtom, radii, tol) and not any_atom_outside:
                    filledAtom.append(new_mol)
                    print("Added molecule:", new_mol)

    return filledAtom

def MoleculesDefiniteBox(numMol, maxattempts, og, box, radii, tol):
    '''Places a defined number of molecules in box'''
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

        if not isOverlapMolecule(new_mol, filledAtom, radii, tol) and not any_atom_outside:
            filledAtom.append(new_mol)
            print("Added molecule:", new_mol)

    return list(filledAtom)

def inBox(xmin, xmax, ymin, ymax, zmin, zmax, box):
    '''Checks if atoms are fully in the box'''
    if (
        xmin >= box.xCoord and xmax <= box.xCoord + box.length and
        ymin >= box.yCoord and ymax <= box.yCoord + box.width and
        zmin >= box.zCoord and zmax <= box.zCoord + box.height
    ):
        return True
    else:
        return False