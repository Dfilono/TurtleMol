import numpy as np
from .makeStruc import calcDistance
from .isOverlap import isOverlapMolecule, isOverlapAtom

def placeMols(shape, og, density, tol, shapeType, radii):
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

    
    mols = []
    for point in gridPoints:
        newPoint = translateMol(og, point)

        if len(og) == 1:
            if not isOverlapAtom(newPoint, mols, radii, tol):
                mols.append(newPoint)
                type = "atom"
        else:
            if not isOverlapMolecule(newPoint, mols, radii, tol):
                mols.append(newPoint)
                type = "molecule"
    
    return mols, type
