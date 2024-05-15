'''Fills an arbitrary mesh with molecules'''

import random
import numpy as np
from .isOverlap import isOverlapAtomKDTree, isOverlapMoleculeKDTree, buildKDTreeMapping
from .makeStruc import makeBase, reCenter, randReorient

def atomsFillMesh(mesh, og, tol, radii, numMol):
    '''Fills mesh with single atoms'''
    filled = []

    if str(numMol).lower() == 'fill':
        numMol = 10000000000000

    # Create KD-tree for filledAtoms
    kdTree, indexToAtom = buildKDTreeMapping(filled, radii)

    # Determine bounds of mesh
    bounds = mesh.bounds
    minBound, maxBound = bounds[0], bounds[1]

    # Determine spacing between molecules
    spacing = tol

    # Generate grid of points
    gridX, gridY, gridZ = np.mgrid[minBound[0]:maxBound[0]:spacing,
                                   minBound[1]:maxBound[1]:spacing,
                                   minBound[2]:maxBound[2]:spacing]

    # Check each point in the grid
    for x in np.nditer(gridX):
        for y in np.nditer(gridY):
            for z in np.nditer(gridZ):
                for atom in og:
                    # Construct atom data
                    atomData = [atom[0], x, y, z]
                    if len(atom) == 5:
                            atomData.append(atom[4])
                    point = [x, y, z]
                    if mesh.isInside(point) and \
                        (kdTree is None or not isOverlapAtomKDTree(atomData, kdTree, indexToAtom, radii, tol)):

                        filled.append(atomData)

                        # Rebuild KDTree with newly added atoms
                        kdTree, indexToAtom = buildKDTreeMapping(filled, radii)
    
    return filled

def moleculesFillMesh(mesh, og, tol, radii, numMol, baseStruc, 
                      randOrient):
    '''Fills mesh with molecules'''
    filled = []

    if str(numMol).lower() == 'fill':
        numMol = 10000000000000

    if baseStruc is not None:
        base = makeBase(baseStruc)
        filled.append(reCenter(base, mesh))

    # Create KD-tree for filledAtoms
    kdTree, indexToAtom = buildKDTreeMapping(filled, radii)

    # Determine bounds of mesh
    bounds = mesh.bounds
    minBound, maxBound = bounds[0], bounds[1]

    # Determine spacing between molecules
    spacing = tol

    # Generate grid of points
    gridX, gridY, gridZ = np.mgrid[minBound[0]:maxBound[0]:spacing,
                                   minBound[1]:maxBound[1]:spacing,
                                   minBound[2]:maxBound[2]:spacing]
    
    # Check each point in the grid
    for i in range(gridX.shape[0]):
        for j in range(gridY.shape[1]):
            for k in range(gridZ.shape[2]):
                x = gridX[i, j, k]
                y = gridY[i, j, k]
                z = gridZ[i, j, k]

                # Check if entire molecule can be placed
                molValid = True
                newMol = []
                for atom in og:
                    # Construct atom data
                    atomType, relX, relY, relZ = atom[:4]
                    atomPoint = [x + relX, y + relY, z + relZ]

                    if not mesh.isInside(atomPoint):
                        molValid = False
                        break
                
                if molValid:
                    for atom in og:
                        atomType, relX, relY, relZ = atom[:4]
                        if len(atom) == 4:
                            atomData = (atomType, float(x + relX), float(y + relY), float(z + relZ))
                        if len(atom) == 5:
                            atomData = (atomType, float(x + relX), float(y + relY), float(z + relZ), atom[4])

                        newMol.append(atomData)
                        
                    if randOrient and len(newMol) == len(og):
                        newMol = randReorient(newMol)
                    if (kdTree is None or not isOverlapMoleculeKDTree(newMol, kdTree, indexToAtom, radii, tol)):
                        filled.append(newMol)

                        # Rebuild KDTree with newly added atoms
                        kdTree, indexToAtom = buildKDTreeMapping(filled, radii)
    return filled

def atomsRandMesh(mesh, og, tol, radii, numMol, maxAttempts):
    '''Randomly places atoms in a mesh'''
    filled = []
    attempts = 0

    # Create KD-tree for filledAtoms
    kdTree, indexToAtom = buildKDTreeMapping(filled, radii)

    # Determine bounds of mesh
    bounds = mesh.bounds
    minBound, maxBound = bounds[0], bounds[1]

    while len(filled) < numMol and attempts <= maxAttempts:
        for atom in og:
            x = random.uniform(minBound[0], maxBound[0])
            y = random.uniform(minBound[1], maxBound[1])
            z = random.uniform(minBound[2], maxBound[2])
            atomType, xRel, yRel, zRel = atom[:4]
            atomPoint = [x + xRel, y + yRel, z + zRel]
            
            if mesh.isInside(atomPoint):
                if len(atom) == 4:
                    atomData = (atomType, atomPoint[0], atomPoint[1], atomPoint[2])

                elif len(atom) == 5:
                    atomData = (atomType, atomPoint[0], atomPoint[1], atomPoint[2], atom[4])

                if (kdTree is None or not isOverlapAtomKDTree(atomData, kdTree, indexToAtom, radii, tol)):
                    filled.append(atomData)

                    # Rebuild KDTree with newly added atoms
                    kdTree, indexToAtom = buildKDTreeMapping(filled, radii)
        attempts += 1
    return filled

def moleculesRandMesh(mesh, og, tol, radii, numMol, baseStruc,
                      randOrient, maxAttempts):
    '''Randomly places molecules in a mesh'''
    filled = []
    attempts = 0

    if baseStruc is not None:
        base = makeBase(baseStruc)
        filled.append(reCenter(base, mesh))

    # Create KD-tree for filledAtoms
    kdTree, indexToAtom = buildKDTreeMapping(filled, radii)

    # Determine bounds of mesh
    bounds = mesh.bounds
    minBound, maxBound = bounds[0], bounds[1]

    while len(filled) < numMol and attempts <= maxAttempts:
        newMol = []

        x = random.uniform(minBound[0], maxBound[0])
        y = random.uniform(minBound[1], maxBound[1])
        z = random.uniform(minBound[2], maxBound[2])

        for atom in og:
            atomType, xRel, yRel, zRel = atom[:4]
            atomPoint = [x + xRel, y + yRel, z + zRel]

            if mesh.isInside(atomPoint):
                if len(atom) == 4:
                    atomData = (atomType, atomPoint[0], atomPoint[1], atomPoint[2])

                elif len(atom) == 5:
                    atomData = (atomType, atomPoint[0], atomPoint[1], atomPoint[2], atom[4])

                newMol.append(atomData)

        if randOrient and len(newMol) == len(og):
            newMol = randReorient(newMol)

        if (kdTree is None or not isOverlapMoleculeKDTree(newMol, kdTree, indexToAtom, radii, tol)):
            if len(newMol) == len(og):
                filled.append(newMol)

                # Rebuild KDTree with newly added atoms
                kdTree, indexToAtom = buildKDTreeMapping(filled, radii)

        attempts += 1
    return list(filled)