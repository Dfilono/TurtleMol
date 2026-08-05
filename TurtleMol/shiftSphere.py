'''Module places molecules in a sphere'''

import random
import scipy.spatial
import numpy as np
from .makeStruc import makeBase, reCenter, Reorient
from .treeManager import PendingKDManager

def atomFillSphere(numShifts, sphere, og, radii, tol, numMol):
    '''Fills sphere with single atoms'''
    filled = []

    if str(numMol).lower() == 'fill':
        numMol = 10000000000000

    # Create KD-tree for filledAtoms
    manager = PendingKDManager(radii, rebuildRate=500)

    for zShifts in range(numShifts):
        for yShifts in range(numShifts):
            for xShifts in range(numShifts):
                newMol = []
                for atom in og:
                    # Calculate the relative position of each atom within the molecule
                    relX = atom[1]
                    relY = atom[2]
                    relZ = atom[3]
                    atomType = atom[0]

                    # Calculate the coordiantes within the sphere
                    x = (sphere.xCoord - sphere.radius) + relX + xShifts*tol
                    y = (sphere.yCoord - sphere.radius) + relY + yShifts*tol
                    z = (sphere.zCoord - sphere.radius) + relZ + zShifts*tol

                    # Adjust for atomic radii
                    atomRadius = radii.get(atomType, 0.0) # Get the radius for the atom type

                    # Check if the new atom fits within the sphere
                    if sphere.containsPoints(x, y, z, atomRadius):
                        if len(atom) == 5:
                            newAtom = (atom[0], x, y, z, atom[4])
                        else:
                            newAtom = (atom[0], x, y, z)

                        newMol.append(newAtom)

                        if (numMol > len(filled) and not manager.overlapsMolecule(newMol, tol)):
                            filled.append(newMol)

                            # Rebuild KDTree with newly added atoms
                            manager.addMolecule(newMol)
                            manager.maybeRebuild() 

                        if len(filled) >= numMol:
                            return filled

    return filled

def atomRandSphere(numMol, maxAttempts, og, sphere, radii, tol):
    '''Places a defined number of atoms in a sphere'''
    filled = []
    attempts = 0

    # Create KD-tree for filledAtoms
    manager = PendingKDManager(radii, rebuildRate=500)

    while len(filled) < numMol and attempts <= maxAttempts:
        newMol = []
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
                if len(atom) == 5:
                    newAtom = (atom[0], newX, newY, newZ, atom[4])
                else:
                    newAtom = (atom[0], newX, newY, newZ)
                
                newMol.append(newAtom)

                if (numMol > len(filled) and not manager.overlapsMolecule(newMol, tol)):
                    filled.append(newMol)

                    # Rebuild KDTree with newly added atoms
                    manager.addMolecule(newMol)
                    manager.maybeRebuild()

                
                if len(filled) >= numMol:
                    return filled

        attempts += 1

    return list(filled)

def moleculeFillSphere(numShifts, sphere, og, radii, tol,
                       baseStruc, randOrient, numMol, rotAngles):
    '''Fills sphere with molecules'''
    filled = []

    if str(numMol).lower() == 'fill':
        numMol = 10000000000000

    # Use pending-buffer KD manager
    manager = PendingKDManager(radii, rebuildRate=500)

    if baseStruc is not None:
        base = makeBase(baseStruc)
        filled.append(reCenter(base, sphere))
        manager.addMolecule(base)
        manager.maybeRebuild(force=True) #ensure the base is comitted

    for zShifts in range(numShifts):
        for yShifts in range(numShifts):
            for xShifts in range(numShifts):
                newMol = []

                for atom in og:
                    # Calculate the relative position of each atom within the molecule
                    relX = atom[1]
                    relY = atom[2]
                    relZ = atom[3]
                    atomType = atom[0]

                    # Calculate the coordiantes within the sphere
                    x = (sphere.xCoord - sphere.radius) + relX + xShifts*tol
                    y = (sphere.yCoord - sphere.radius) + relY + yShifts*tol
                    z = (sphere.zCoord - sphere.radius) + relZ + zShifts*tol

                    # Adjust for atomic radii
                    atomRadius = radii.get(atomType, 0.0) # Get the radius for the atom type

                    # Check if the new atom fits within the sphere
                    if sphere.containsPoints(x, y, z, atomRadius):
                        if len(atom) == 5:
                            newAtom = (atom[0], x, y, z, atom[4])
                        else:
                            newAtom = (atom[0], x, y, z)
                        newMol.append(newAtom)
                    else:
                        break # If any atom doesn't fit, discard the whol molecule

                if randOrient and len(newMol) == len(og):
                    newMol = Reorient(newMol, randRotate=True)

                if not np.allclose(np.array(rotAngles), np.array([0, 0, 0])) and len(newMol) == len(og):
                    newMol = Reorient(newMol, angles=rotAngles)

                if (numMol > len(filled) and not manager.overlapsMolecule(newMol, tol)):
                    if len(newMol) == len(og):
                        filled.append(newMol)

                        # Rebuild KDTree with newly added atoms
                        manager.addMolecule(newMol)
                        manager.maybeRebuild()

                
                if len(filled) >= numMol:
                    return filled
                
    return filled

def moleculeRandSphere(numMol, maxAttempts, og, sphere, radii, tol,
                           baseStruc, randOrient, rotAngles):
    '''Places a defined number of molecules in a sphere'''
    filled = []
    attempts = 0

    # Use pending-buffer KD manager
    manager = PendingKDManager(radii, rebuildRate=500)

    if baseStruc is not None:
        base = makeBase(baseStruc)
        filled.append(reCenter(base, sphere))
        manager.addMolecule(base)
        manager.maybeRebuild(force=True) #ensure the base is comitted

    while len(filled) < numMol and attempts <= maxAttempts:
        newMol = []

        shiftX = 0
        shiftY = 0
        shiftZ = 0

        shiftX = random.uniform((sphere.xCoord - sphere.radius),
                                (sphere.xCoord + sphere.radius))
        shiftY = random.uniform((sphere.yCoord - sphere.radius),
                                (sphere.yCoord + sphere.radius))
        shiftZ = random.uniform((sphere.zCoord - sphere.radius),
                                    (sphere.zCoord + sphere.radius))
        for atom in og:

            # Calculate the shift for each tile and new point
            newX = atom[1] + shiftX
            newY = atom[2] + shiftY
            newZ = atom[3] + shiftZ
            atomType = atom[0]

            # Adjust for atomic radii
            atomRadius = radii.get(atomType, 0.0) # Get the radius for the atom type
            # Check if the new atom fits within the sphere
            if sphere.containsPoints(newX, newY, newZ, atomRadius):
                if len(atom) == 5:
                    newAtom = (atom[0], newX, newY, newZ, atom[4])
                else:
                    newAtom = (atom[0], newX, newY, newZ)
                newMol.append(newAtom)

                if randOrient and len(newMol) == len(og):
                    newMol = Reorient(newMol, randRotate=True)

                if not np.allclose(np.array(rotAngles), np.array([0, 0, 0])) and len(newMol) == len(og):
                    newMol = Reorient(newMol, angles=rotAngles)

                if (numMol > len(filled) and not manager.overlapsMolecule(newMol, tol)):
                    if len(newMol) == len(og):
                        filled.append(newMol)

                        # Rebuild KDTree with newly added atoms
                        manager.addMolecule(newMol)
                        manager.maybeRebuild()

        
        if len(filled) >= numMol:
            return filled

        attempts += 1

    return list(filled)
