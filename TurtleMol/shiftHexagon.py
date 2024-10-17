import numpy as np
from .shiftBox import inBox

def hexagonUnitCellBox(shape, dims, cellDims, og, radii):
    '''Duplicates hexagonal unit cells to fill a given box'''

    # Find the atom types in the tile
    atomNames = {atom[0] for atom in og}
    totalRadius = sum(radii[str.capitalize(name)] for name in atomNames)

    # Calculate how many times to duplicate the unit cell in a given dimension
    dupeCount = [int(dims[i] / cellDims[i]) for i in range(3)]
    cellParams = f'CRYST1{dupeCount[0]*cellDims[0]:9.3f}{dupeCount[1]*cellDims[1]:9.3f}{dupeCount[2]*cellDims[2]:9.3f}  90.00  90.00  120.0 P 1         1'

    # Define the basis vectors for hexagonal lattice
    a1 = np.array([cellDims[0], 0, 0])
    a2 = np.array([cellDims[1] * np.cos(np.radians(120)), cellDims[1] * np.sin(np.radians(120)), 0])
    a3 = np.array([0, 0, cellDims[2]])

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

def hexagonUnitCellSphere(shape, cellDims, og, radii):
    '''Duplicates unit cells to fill a given sphere'''

    # Box dimensions that completely contain the sphere
    boxDim = [2 * shape.radius] * 3
    dupeCount = [int(boxDim[i] / cellDims[i]) for i in range(3)]
    cellParams = f'CRYST1{dupeCount[0]*cellDims[0]:9.3f}{dupeCount[1]*cellDims[1]:9.3f}{dupeCount[2]*cellDims[2]:9.3f}  90.00  90.00  120.0 P 1         1'

    # Define the basis vectors for hexagonal lattice
    a1 = np.array([cellDims[0], 0, 0])
    a2 = np.array([cellDims[1] * np.cos(np.radians(120)), cellDims[1] * np.sin(np.radians(120)), 0])
    a3 = np.array([0, 0, cellDims[2]])


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

def hexagonUnitCellMesh(shape, cellDims, og, radius):
    '''Duplicates unit cells to fill a given mesh'''

    # Find the atom types in the tile
    atomNames = {atom[0] for atom in og}
    totalRadius = sum(radius[str.capitalize(name)] for name in atomNames)

    # Box dimensions that completely contain the mesh
    maxBound, minBound = shape.bounds[1], shape.bounds[0]
    boxDim = maxBound - minBound
    dupeCount = [int(boxDim[i] / cellDims[i]) for i in range(3)]
    
    cellParams = f'CRYST1{dupeCount[0]*cellDims[0]:9.3f}{dupeCount[1]*cellDims[1]:9.3f}{dupeCount[2]*cellDims[2]:9.3f}  90.00  90.00  120.0 P 1         1'

    # Define the basis vectors for hexagonal lattice
    a1 = np.array([cellDims[0], 0, 0])
    a2 = np.array([cellDims[1] * np.cos(np.radians(120)), cellDims[1] * np.sin(np.radians(120)), 0])
    a3 = np.array([0, 0, cellDims[2]])

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