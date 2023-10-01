'''Module places molecules in a box'''

import random
from isOverlap import isOverlapAtom, isOverlapMolecule

def atomsFillBox(x, y, z, tol, og, box, radii):
    '''Fills box with single atoms'''
    filledAtom = []
    filledPositions = set()

    for zShift in range(z):
        for yShift in range(y):
            for xShift in range(x):
                for atom in og:
                    # Calculate the shift for each tile and new point
                    newX = box.xCoord + atom[1] + xShift
                    newY = box.yCoord + atom[2] + yShift
                    newZ = box.zCoord + atom[3] + zShift

                    # Adjust for atomic radiss
                    newXMin = newX - radii[atom[0]]
                    newYMin = newY - radii[atom[0]]
                    newZMin = newZ - radii[atom[0]]
                    newXMax = newX + radii[atom[0]]
                    newYMax = newY + radii[atom[0]]
                    newZMax = newZ + radii[atom[0]]

                    # Check if the new atom fits within the box
                    if inBox(newXMin, newXMax, newYMin, newYMax, newZMin, newZMax, box):
                        newAtom = (atom[0], newX, newY, newZ)

                        if not isOverlapAtom(newAtom, filledAtom, radii, tol):
                            filledAtom.append(newAtom)
                            filledPositions.add(newAtom)

    return filledAtom

def atomsDefiniteBox(numMol, maxAttempts, og, box, radii, tol):
    '''Places a defined number of single atoms in box'''
    filledAtom = []
    filledPositions = set()

    attempts = 0

    while len(filledAtom) < numMol and attempts < maxAttempts:
        attempts += 1
        for atom in og:

            # Calculate the shift for each tile and new point
            newX = atom[1] + random.uniform(0, box.length)
            newY = atom[2] + random.uniform(0, box.width)
            newZ = atom[3] + random.uniform(0, box.height)

            # Adjust for atomic radiss
            newXMin = newX - radii[atom[0]]
            newYMin = newY - radii[atom[0]]
            newZMin = newZ - radii[atom[0]]
            newXMax = newX + radii[atom[0]]
            newYMax = newY + radii[atom[0]]
            newZMax = newZ + radii[atom[0]]

            # Check if the new atom fits within the box
            if  inBox(newXMin, newXMax, newYMin, newYMax, newZMin, newZMax, box):
                newAtom = (atom[0], newX, newY, newZ)

                if not isOverlapAtom(newAtom, filledAtom, radii, tol):
                    filledAtom.append(newAtom)
                    filledPositions.add(newAtom)

    return list(filledAtom)

def moleculesFillBox(x, y, z, tol, og, box, radii):
    '''Fills box with molecules'''
    filledAtom = []

    for zShift in range(z):
        for yShift in range(y):
            for xShift in range(x):
                newMol = []
                anyAtomOutside = False

                for i, atom in enumerate(og):
                    #print(atom)
                    # Calculate the shift for each tile and new point
                    newX = box.xCoord + atom[1] + xShift
                    newY = box.yCoord + atom[2] + yShift
                    newZ = box.zCoord + atom[3] + zShift

                    # Adjust for atomic radiss
                    newXMin = newX - radii[atom[0]]
                    newYMin = newY - radii[atom[0]]
                    newZMin = newZ - radii[atom[0]]
                    newXMax = newX + radii[atom[0]]
                    newYMax = newY + radii[atom[0]]
                    newZMax = newZ + radii[atom[0]]

                    # Check if the new atom fits within the box
                    if  inBox(newXMin, newXMax, newYMin, newYMax, newZMin,
                              newZMax, box):
                        newAtom = (atom[0], newX, newY, newZ)
                        newMol.append(newAtom)

                    else:
                        anyAtomOutside = True

                if (not isOverlapMolecule(newMol, filledAtom, radii, tol) and
                    not anyAtomOutside):

                    filledAtom.append(newMol)
                    print("Added molecule:", newMol)

    return filledAtom

def moleculesDefiniteBox(numMol, maxAttempts, og, box, radii, tol):
    '''Places a defined number of molecules in box'''
    filledAtom = []

    attempts = 0

    while len(filledAtom) < numMol and attempts < maxAttempts:
        attempts += 1
        newMol = []
        anyAtomOutside = False

        xShift = random.uniform(0, box.length)
        yShift = random.uniform(0, box.width)
        zShift = random.uniform(0, box.height)

        for atom in og:

            # Calculate the shift for each tile and new point
            newX = atom[1] + xShift
            newY = atom[2] + yShift
            newZ = atom[3] + zShift

            # Adjust for atomic radiss
            newXMin = newX - radii[atom[0]]
            newYMin = newY - radii[atom[0]]
            newZMin = newZ - radii[atom[0]]
            newXMax = newX + radii[atom[0]]
            newYMax = newY + radii[atom[0]]
            newZMax = newZ + radii[atom[0]]

            # Check if the new atom fits within the box
            if  inBox(newXMin, newXMax, newYMin, newYMax, newZMin, newZMax, box):
                newAtom = (atom[0], newX, newY, newZ)
                newMol.append(newAtom)

            else:
                anyAtomOutside = True

        if not isOverlapMolecule(newMol, filledAtom, radii, tol) and not anyAtomOutside:
            filledAtom.append(newMol)
            print("Added molecule:", newMol)

    return list(filledAtom)

def inBox(xMin, xMax, yMin, yMax, zMin, zMax, box):
    '''Checks if atoms are fully in the box'''
    if (
        xMin >= box.xCoord and xMax <= box.xCoord + box.length and
        yMin >= box.yCoord and yMax <= box.yCoord + box.width and
        zMin >= box.zCoord and zMax <= box.zCoord + box.height
    ):
        return True
    return False