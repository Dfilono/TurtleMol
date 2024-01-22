'''This module is for replicating a unit cell if it is defined in a given pdb'''

import numpy as np
import trimesh
from .shiftBox import inBox

# Box
def unitCellBox(shape, dims, cellDims, og, radii):
    '''Duplicates unit cells to fill a given box'''
    # Calculate how many times to duplicate the unit cell in a given dimension
    dupeCount = [int(dims[i] / cellDims[i]) for i in range(3)]

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
                    newX = atom[1] + disp
                    newY = atom[2] + disp
                    newZ = atom[3] + disp
                    atomType = atom[0]

                    # Adjust for atomic radiss
                    newXMin = newX - radii[atom[0]]
                    newYMin = newY - radii[atom[0]]
                    newZMin = newZ - radii[atom[0]]
                    newXMax = newX + radii[atom[0]]
                    newYMax = newY + radii[atom[0]]
                    newZMax = newZ + radii[atom[0]]

                    # Check if the new atom fits within the box
                    if  inBox(newXMin, newXMax, newYMin, newYMax, newZMin,
                              newZMax, shape):
                        if len(atom) == 5:
                            newAtom = (atom[0], newX, newY, newZ, atom[4])
                        else:
                            newAtom = (atom[0], newX, newY, newZ)
                        currentCell.append(newAtom)

                filled.append(currentCell)
    return filled, "molecule"

def unitCellSphere(shape, cellDims, og, radii):
    '''Duplicates unit cells to fill a given sphere'''
    # Box dimensions that completely contain the sphere
    boxDim = [2 * shape.radius] * 3
    dupeCount = [int(boxDim[i] / cellDims[i]) for i in range(3)]

    filled = []

    for dx in range(dupeCount[0]):
        for dy in range(dupeCount[1]):
            for dz in range(dupeCount[2]):
                disp = np.array([dx * cellDims[0], dy * cellDims[1], dz * cellDims[2]])
                currentCell = []

                for atom in og:
                    newX = atom[1] + disp
                    newY = atom[2] + disp
                    newZ = atom[3] + disp
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
    return filled, "molecule"

def unitCellMesh(shape, cellDims, og):
    '''Duplicates unit cells to fill a given mesh'''
    # Box dimensions that completely contain the mesh
    maxBound, minBound = shape.bounds[0], shape.bounds[1]
    boxDim = maxBound - minBound
    dupeCount = [int(boxDim[i] / cellDims[i]) for i in range(3)]

    filled = []

    for dx in range(dupeCount[0]):
        for dy in range(dupeCount[1]):
            for dz in range(dupeCount[2]):
                disp = np.array([dx * cellDims[0], dy * cellDims[1], dz * cellDims[2]])
                currentCell = []

                for atom in og:
                    newX = atom[1] + disp
                    newY = atom[2] + disp
                    newZ = atom[3] + disp
                    atomType = atom[0]

                    if shape.isInside([newX, newY, newZ]):
                        if len(atom) == 5:
                            newAtom = (atom[0], newX, newY, newZ, atom[4])
                        else:
                            newAtom = (atom[0], newX, newY, newZ)
                        currentCell.append(newAtom)
                filled.append(currentCell)
    return filled, "molecule"