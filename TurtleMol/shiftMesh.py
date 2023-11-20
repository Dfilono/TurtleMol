'''Fills an arbitrary mesh with molecules'''

import numpy as np
import random
from .isOverlap import isOverlapAtom, isOverlapMolecule
from .makeStruc import makeBase, reCenter, randReorient

def atomsFillMesh(mesh, og, tol, radii, numMol):
    '''Fills mesh with single atoms'''
    filled = []

    if str(numMol).lower() == 'fill':
        numMol = 10000000000000

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
                        not isOverlapAtom(atomData, filled, radii, tol):

                        filled.append(atomData)
    
    return filled

def moleculesFillMesh(mesh, og, tol, radii, numMol, baseStruc):
    '''Fills mesh with molecules'''
    filled = []

    if str(numMol).lower() == 'fill':
        numMol = 10000000000000

    if baseStruc is not None:
        base = makeBase(baseStruc)
        filled.append(reCenter(base, ))

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
                        not isOverlapMolecule(atomData, filled, radii, tol):
                        
                        filled.append(atomData)
    
    return filled
