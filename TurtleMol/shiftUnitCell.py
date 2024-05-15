'''This module is for replicating a unit cell if it is defined in a given pdb'''

import numpy as np
import trimesh
from .shiftBox import inBox

# Box
def unitCellBox(shape, dims, cellDims, og, radii):
    '''Duplicates unit cells to fill a given box'''

    # Find the atom types in the tile
    atomNames = {atom[0] for atom in og}
    totalRadius = sum(radii[str.capitalize(name)] for name in atomNames)

    # Calculate how many times to duplicate the unit cell in a given dimension
    dupeCount = [int(dims[i] / cellDims[i]) for i in range(3)]
    cellParams = f'CRYST1    {dupeCount[0]*cellDims[0]: .3f}    {dupeCount[1]*cellDims[1]: .3f}    {dupeCount[2]*cellDims[2]: .3f}  90.00  90.00  90.00 P1          1'

    filled = []

    for dx in range(dupeCount[0]):
        for dy in range(dupeCount[1]):
            for dz in range(dupeCount[2]):
                # Calculate the displacement for this duplication
                disp = np.array([dx * cellDims[0], 
                                 dy * cellDims[1], 
                                 dz * cellDims[2]])
                
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

def unitCellSphere(shape, cellDims, og, radii):
    '''Duplicates unit cells to fill a given sphere'''

    # Box dimensions that completely contain the sphere
    boxDim = [2 * shape.radius] * 3
    dupeCount = [int(boxDim[i] / cellDims[i]) for i in range(3)]
    cellParams = f'CRYST1    {dupeCount[0]*cellDims[0]: .3f}    {dupeCount[1]*cellDims[1]: .3f}    {dupeCount[2]*cellDims[2]: .3f}  90.00  90.00  90.00 P1          1'

    filled = []

    for dx in range(dupeCount[0]):
        for dy in range(dupeCount[1]):
            for dz in range(dupeCount[2]):
                disp = np.array([dx * cellDims[0], 
                                 dy * cellDims[1], 
                                 dz * cellDims[2]])
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

def unitCellMesh(shape, cellDims, og, radius):
    '''Duplicates unit cells to fill a given mesh'''

    # Find the atom types in the tile
    atomNames = {atom[0] for atom in og}
    totalRadius = sum(radius[str.capitalize(name)] for name in atomNames)

    # Box dimensions that completely contain the mesh
    maxBound, minBound = shape.bounds[1], shape.bounds[0]
    boxDim = maxBound - minBound
    dupeCount = [int(boxDim[i] / cellDims[i]) for i in range(3)]
    
    cellParams = f'CRYST1    {dupeCount[0]*cellDims[0]: .3f}    {dupeCount[1]*cellDims[1]: .3f}    {dupeCount[2]*cellDims[2]: .3f}  90.00  90.00  90.00 P1          1'

    filled = []

    for dx in range(dupeCount[0]):
        for dy in range(dupeCount[1]):
            for dz in range(dupeCount[2]):
                disp = np.array([dx * cellDims[0], 
                                 dy * cellDims[1], 
                                 dz * cellDims[2]]) + minBound
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