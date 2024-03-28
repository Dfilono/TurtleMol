'''Places molecules in a given volume based on the provided denisty'''

import numpy as np
import scipy.spatial
from .makeStruc import calcDistance, randReorient
from .isOverlap import isOverlapMoleculeKDTree, isOverlapAtomKDTree

def placeMols(shape, og, density, tol, shapeType, radii, randOrient):
    '''Place molecules in the grid defined by the density'''
    gridPoints = calcDistance(shape, og, density, shapeType)

    def calcCentroid(coords):
        """Calculate the geometric centroid of a set of coordinates."""
        return np.mean(coords, axis=0)

    # Tanslate molecules such that their centroid is at the center
    # of the target point
    def translateMol(coords, target):
        '''Translate molecules to a new grid point'''
        atomLabels = [atom[0] for atom in coords]
        pdbInfo = [atom[4] if len(atom) == 5 else None for atom in coords]
        xyzCoords = np.array([atom[1:4] for atom in coords])

        centroid = calcCentroid(xyzCoords)
        transVec = np.array(target) - centroid
        translatedXYZ = xyzCoords + transVec

        transMol = [[label] + coord.tolist() + ([info] if info is not None else [])
                    for label, coord, info in zip(atomLabels, translatedXYZ, pdbInfo)]
        
        return transMol

    strucType = 'molecule'
    mols = []

    # Create KD-tree for filledAtoms
    kdTree = scipy.spatial.cKDTree(np.array(mols))

    for point in gridPoints:
        newPoint = translateMol(og, point)

        if randOrient:
            newPoint = randReorient(newPoint)

        if len(og) == 1:
            if not isOverlapAtomKDTree(newPoint, kdTree, radii, tol):
                mols.append(newPoint)

                # Rebuild KDTree with newly added atoms
                kdTree = scipy.spatial.cKDTree(np.array([atom[1:4] for atom in mols]))

                strucType = "atom"
        else:
            if not isOverlapMoleculeKDTree(newPoint, kdTree, radii, tol):
                mols.append(newPoint)

                # Rebuild KDTree with newly added atoms
                kdTree = scipy.spatial.cKDTree(np.array([atom[1:4] for atom in mols]))

                strucType = "molecule"
    
    return mols, strucType
