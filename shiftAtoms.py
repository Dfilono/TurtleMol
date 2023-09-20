import random
import math

def shiftAtomsFill(x, y, z, tol, og, box, radii):
    filledAtom = []
    filledPositions = set()
    for zshift in range(z):
        for yshift in range(y):
            for xshift in range(x):
                for atom in og:
                    # Calculate the shift for each tile and new point
                    new_x = atom[1] + xshift * tol
                    new_y = atom[2] + yshift * tol
                    new_z = atom[3] + zshift * tol

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

                if not is_overlap(new_atom, filledAtom, radii):
                    filledAtom.append(new_atom)
                    filledPositions.add(new_atom)

    return filledAtom

def shiftAtomsDefinite(numMol, maxattempts, og, box, radii): # NOTE only works by randomly assigning atoms
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
            if  (
                new_x_min >= box.x and new_x_max <= box.x + box.length and
                new_y_min >= box.y and new_y_max <= box.y + box.width and
                new_z_min >= box.z and new_z_max <= box.z + box.height
            ):
                new_atom = (atom[0], new_x, new_y, new_z)

                if not is_overlap(new_atom, filledAtom, radii):
                    filledAtom.append(new_atom)
                    filledPositions.add(new_atom)


    return list(filledAtom)

def is_overlap(new_atom, filled_atoms, radii):
    for atom in filled_atoms:
        distance = math.sqrt(
            (new_atom[1] - atom[1])**2 +
            (new_atom[2] - atom[2])**2 +
            (new_atom[3] - atom[3])**2
        )

        if distance < (radii[new_atom[0]] + radii[atom[0]]):
            #print(f'Fail: dist = {distance}, Radius = {radii[new_atom[0]] + radii[atom[0]]}')
            return True
        #print(f'Pass: dist = {distance}, Radius = {radii[new_atom[0]] + radii[atom[0]]}')
    return False