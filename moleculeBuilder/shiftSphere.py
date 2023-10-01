'''Module places molecules in a sphere'''

import random
from isOverlap import isOverlapAtom, isOverlapMolecule

def atomFillSphere(numShifts, sphere, og, radii, tol):
    '''Fills sphere with single atoms'''
    filled = []

    for zShifts in range(numShifts):
        for yShifts in range(numShifts):
            for xShifts in range(numShifts):
                for atom in og:
                    # Calculate the relative position of each atom within the molecule
                    relX = atom[1]
                    relY = atom[2]
                    relZ = atom[3]
                    atomType = atom[0]

                    # Calculate the coordiantes within the sphere
                    x = (sphere.xCoord - sphere.radius) + relX + xShifts
                    y = (sphere.yCoord - sphere.radius) + relY + yShifts
                    z = (sphere.zCoord - sphere.radius) + relZ + zShifts

                    # Adjust for atomic radii
                    atomRadius = radii.get(atomType, 0.0) # Get the radius for the atom type

                    # Check if the new atom fits within the sphere
                    if sphere.containsPoints(x, y, z, atomRadius):
                        newAtom = (atomType, x, y, z)

                        if not isOverlapAtom(newAtom, filled, radii, tol):
                            filled.append(newAtom)

    return filled

def atomDefiniteSphere(numMol, maxAttempts, og, sphere, radii, tol):
    '''Places a defined number of atoms in a sphere'''
    filled = []
    attempts = 0

    while len(filled) < numMol and attempts <= maxAttempts:
        for atom in og:

            # Calculate the shift for each tile and new point
            newX = atom[1] + random.uniform((sphere.xCoord - sphere.radius),
                                            (sphere.xCoord + sphere.radius))
            newY = atom[2] + random.uniform((sphere.yCoord - sphere.radius),
                                            (sphere.yCoord + sphere.radius))
            newZ = atom[3] + random.uniform((sphere.zCoord - sphere.radius),
                                            (sphere.zCoord + sphere.radius))
            atomType = atom[0]

            # Adjust for atomic radii
            atomRadius = radii.get(atomType, 0.0) # Get the radius for the atom type

            # Check if the new atom fits within the sphere
            if sphere.containsPoints(newX, newY, newZ, atomRadius):
                newAtom = (atomType, newX, newY, newZ)

                if not isOverlapAtom(newAtom, filled, radii, tol):
                    filled.append(newAtom)

        attempts += 1

    return list(filled)

def moleculeFillSphere(numShifts, sphere, og, radii, tol):
    '''Fills sphere with molecules'''
    filled = []

    for zShifts in range(numShifts):
        for yShifts in range(numShifts):
            for xShifts in range(numShifts):
                newMol = []

                for i, atom in enumerate(og):
                    # Calculate the relative position of each atom within the molecule
                    relX = atom[1]
                    relY = atom[2]
                    relZ = atom[3]
                    atomType = atom[0]

                    # Calculate the coordiantes within the sphere
                    x = (sphere.xCoord - sphere.radius) + relX + xShifts
                    y = (sphere.yCoord - sphere.radius) + relY + yShifts
                    z = (sphere.zCoord - sphere.radius) + relZ + zShifts

                    # Adjust for atomic radii
                    atomRadius = radii.get(atomType, 0.0) # Get the radius for the atom type

                    # Check if the new atom fits within the sphere
                    if sphere.containsPoints(x, y, z, atomRadius):
                        newAtom = (atomType, x, y, z)
                        newMol.append(newAtom)
                    else:
                        break # If any atom doesn't fit, discard the whol molecule

                if not isOverlapMolecule(newMol, filled, radii, tol):
                    if len(newMol) == len(og):
                        filled.append(newMol)

    return filled

def moleculeDefiniteSphere(numMol, maxAttempts, og, sphere, radii, tol):
    '''Places a defined number of molecules in a sphere'''
    filled = []
    attempts = 0

    while len(filled) < numMol and attempts <= maxAttempts:
        newMol = []
        shiftx = random.uniform((sphere.xCoord - sphere.radius),
                                (sphere.xCoord + sphere.radius))
        shifty = random.uniform((sphere.yCoord - sphere.radius),
                                (sphere.yCoord + sphere.radius))
        shiftz = random.uniform((sphere.zCoord - sphere.radius),
                                (sphere.zCoord + sphere.radius))

        for atom in og:

            # Calculate the shift for each tile and new point
            newX = atom[1] + shiftx
            newY = atom[2] + shifty
            newZ = atom[3] + shiftz
            atomType = atom[0]

            # Adjust for atomic radii
            atomRadius = radii.get(atomType, 0.0) # Get the radius for the atom type
            # Check if the new atom fits within the sphere
            if sphere.containsPoints(newX, newY, newZ, atomRadius):
                newAtom = (atomType, newX, newY, newZ)
                newMol.append(newAtom)
            else:
                break # If any atom doesn't fit, discard the whol molecule

        if not isOverlapMolecule(newMol, filled, radii, tol):
            if len(newMol) == len(og):
                filled.append(newMol)

        attempts += 1

    return list(filled)
