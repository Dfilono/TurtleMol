'''This module is for replicating a unit cell if it is defined in a given pdb'''

import numpy as np
import trimesh
from .shiftBox import inBox
from .makeStruc import calcLatticeVectors, rotateUnitCell

# Box
def unitCellBox(shape, dims, cellDims, cellAngles, og, radii, rotAngles=np.array([0, 0, 0]).all()):
    '''Duplicates unit cells to fill a given box'''

    # Find the atom types in the tile
    atomNames = {atom[0] for atom in og}
    totalRadius = sum(radii[str.capitalize(name)] for name in atomNames)

    if np.array(rotAngles).all() != np.array([0, 0, 0]).all():
        latticeVec = calcLatticeVectors(cellDims[0], cellDims[1], cellDims[2], cellAngles[0], cellAngles[1], cellAngles[2])
        og, cellInfo = rotateUnitCell(latticeVec, og, rotAngles)
        cellDims = [cellInfo['a'], cellInfo['b'], cellInfo['c']]
        cellAngles = [cellInfo['alpha'], cellInfo['beta'], cellInfo['gamma']]

    alpha = np.radians(cellAngles[0])
    beta = np.radians(cellAngles[1])
    gamma = np.radians(cellAngles[2])

    # Define the basis vectors
    a1 = np.array([cellDims[0], 0, 0])
    a2 = np.array([cellDims[1] * np.cos(gamma), cellDims[1] * np.sin(gamma), 0])
    a3 = np.array([
        cellDims[2] * np.cos(beta), 
        cellDims[2] * (np.cos(alpha) - np.cos(beta) * np.cos(gamma)) / np.sin(gamma),
        cellDims[2] * np.sqrt(1 - np.cos(beta)**2 - ((np.cos(alpha) - np.cos(beta) * np.cos(gamma)) / np.sin(gamma))**2)
    ])

    # Calculate how many times to duplicate the unit cell in a given dimension
    dupeCount = [int(dims[i] / cellDims[i]) for i in range(3)]
    cellParams = f'CRYST1{dupeCount[0]*cellDims[0]:9.3f}{dupeCount[1]*cellDims[1]:9.3f}{dupeCount[2]*cellDims[2]:9.3f}  {cellAngles[0]:0.2f}  {cellAngles[1]:0.2f}  {cellAngles[2]:0.2f} P 1         1'

    filled = []

    for dx in range(dupeCount[0]):
        for dy in range(dupeCount[1]):
            for dz in range(dupeCount[2]):
                # Calculate the displacement for this duplication
                disp = dx*a1 + dy*a2 + dz*a3
                
                currentCell = []
                
                for atom in og:
                    newX = shape.xCoord + atom[1] + disp[0]
                    newY = shape.yCoord + atom[2] + disp[1]
                    newZ = shape.zCoord + atom[3] + disp[2]
                    atomType = atom[0]

                    atomRadius = radii.get(atomType, 0.0)

                    # Adjust for atomic radiss
                    newXMin = newX - atomRadius
                    newYMin = newY - atomRadius
                    newZMin = newZ - atomRadius
                    newXMax = newX + atomRadius
                    newYMax = newY + atomRadius
                    newZMax = newZ + atomRadius

                    # Check if the new atom fits within the box
                    if  inBox(newXMin, newXMax, newYMin, newYMax, newZMin,
                              newZMax, shape):
                        if len(atom) == 5:
                            newAtom = (atom[0], newX, newY, newZ, atom[4])
                        else:
                            newAtom = (atom[0], newX, newY, newZ)
                        currentCell.append(newAtom)

                filled.append(currentCell)
    return filled, "molecule", cellParams

def unitCellSphere(shape, cellDims, cellAngles, og, radii, rotAngles=np.array([0, 0, 0]).all()):
    '''Duplicates unit cells to fill a given sphere'''

    # Box dimensions that completely contain the sphere
    boxDim = [2 * shape.radius] * 3

    if np.array(rotAngles).all() != np.array([0, 0, 0]).all():
        latticeVec = calcLatticeVectors(cellDims[0], cellDims[1], cellDims[2], cellAngles[0], cellAngles[1], cellAngles[2])
        og, cellInfo = rotateUnitCell(latticeVec, og, rotAngles)
        cellDims = [cellInfo['a'], cellInfo['b'], cellInfo['c']]
        cellAngles = [cellInfo['alpha'], cellInfo['beta'], cellInfo['gamma']]

    alpha = np.radians(cellAngles[0])
    beta = np.radians(cellAngles[1])
    gamma = np.radians(cellAngles[2])

    dupeCount = [int(boxDim[i] / cellDims[i]) for i in range(3)]
    cellParams = f'CRYST1    {dupeCount[0]*cellDims[0]: .3f}    {dupeCount[1]*cellDims[1]: .3f}    {dupeCount[2]*cellDims[2]: .3f}  {cellAngles[0]:0.2f}  {cellAngles[1]:0.2f}  {cellAngles[2]:0.2f} P1          1'

    # Define the basis vectors
    a1 = np.array([cellDims[0], 0, 0])
    a2 = np.array([cellDims[1] * np.cos(gamma), cellDims[1] * np.sin(gamma), 0])
    a3 = np.array([
        cellDims[2] * np.cos(beta), 
        cellDims[2] * (np.cos(alpha) - np.cos(beta) * np.cos(gamma)) / np.sin(gamma),
        cellDims[2] * np.sqrt(1 - np.cos(beta)**2 - ((np.cos(alpha) - np.cos(beta) * np.cos(gamma)) / np.sin(gamma))**2)
    ])

    filled = []

    for dx in range(dupeCount[0]):
        for dy in range(dupeCount[1]):
            for dz in range(dupeCount[2]):
                disp = dx*a1 + dy*a2 + dz*a3
                currentCell = []

                for atom in og:
                    newX = (shape.xCoord - shape.radius)+ atom[1] + disp[0]
                    newY = (shape.yCoord - shape.radius) + atom[2] + disp[1]
                    newZ = (shape.zCoord - shape.radius) + atom[3] + disp[2]
                    atomType = atom[0]

                    atomRadius = radii.get(atomType, 0.0)

                    # Check if the new atom fits within the sphere
                    if shape.containsPoints(newX, newY, newZ, atomRadius):
                        if len(atom) == 5:
                            newAtom = (atom[0], newX, newY, newZ, atom[4])
                        else:
                            newAtom = (atom[0], newX, newY, newZ)
                        currentCell.append(newAtom)
                filled.append(currentCell)
    return filled, "molecule", cellParams

def unitCellMesh(shape, cellDims, cellAngles, og, radius, rotAngles=np.array([0, 0, 0]).all()):
    '''Duplicates unit cells to fill a given mesh'''

    # Find the atom types in the tile
    atomNames = {atom[0] for atom in og}
    totalRadius = sum(radius[str.capitalize(name)] for name in atomNames)

    # Box dimensions that completely contain the mesh
    maxBound, minBound = shape.bounds[1], shape.bounds[0]
    boxDim = maxBound - minBound

    if np.array(rotAngles).all() != np.array([0, 0, 0]).all():
        latticeVec = calcLatticeVectors(cellDims[0], cellDims[1], cellDims[2], cellAngles[0], cellAngles[1], cellAngles[2])
        og, cellInfo = rotateUnitCell(latticeVec, og, rotAngles)
        cellDims = [cellInfo['a'], cellInfo['b'], cellInfo['c']]
        cellAngles = [cellInfo['alpha'], cellInfo['beta'], cellInfo['gamma']]

    alpha = np.radians(cellAngles[0])
    beta = np.radians(cellAngles[1])
    gamma = np.radians(cellAngles[2])

    dupeCount = [int(boxDim[i] / cellDims[i]) for i in range(3)]
    
    cellParams = f'CRYST1    {dupeCount[0]*cellDims[0]: .3f}    {dupeCount[1]*cellDims[1]: .3f}    {dupeCount[2]*cellDims[2]: .3f}  {cellAngles[0]:0.2f}  {cellAngles[1]:0.2f}  {cellAngles[2]:0.2f} P1          1'

    # Define the basis vectors
    a1 = np.array([cellDims[0], 0, 0])
    a2 = np.array([cellDims[1] * np.cos(gamma), cellDims[1] * np.sin(gamma), 0])
    a3 = np.array([
        cellDims[2] * np.cos(beta), 
        cellDims[2] * (np.cos(alpha) - np.cos(beta) * np.cos(gamma)) / np.sin(gamma),
        cellDims[2] * np.sqrt(1 - np.cos(beta)**2 - ((np.cos(alpha) - np.cos(beta) * np.cos(gamma)) / np.sin(gamma))**2)
    ])

    filled = []

    for dx in range(dupeCount[0]):
        for dy in range(dupeCount[1]):
            for dz in range(dupeCount[2]):
                disp = disp = dx*a1 + dy*a2 + dz*a3 + minBound
                currentCell = []

                for atom in og:
                    newX = shape.origin()[0] + atom[1] + disp[0]
                    newY = shape.origin()[1] + atom[2] + disp[1]
                    newZ = shape.origin()[2] + atom[3] + disp[2]
                    atomType = atom[0]

                    if shape.isInside([newX, newY, newZ]):
                        if len(atom) == 5:
                            newAtom = (atom[0], newX, newY, newZ, atom[4])
                        else:
                            newAtom = (atom[0], newX, newY, newZ)
                        currentCell.append(newAtom)
                filled.append(currentCell)
    return filled, "molecule", cellParams
