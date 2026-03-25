'''This module is for replicating a unit cell if it is defined in a given pdb'''

import numpy as np
import math
import trimesh
from .shiftBox import inBox
from .makeStruc import calcLatticeVectors, rotateUnitCell

# Box

def unitCellBox(shape, dims, cellDims, cellAngles, og, radii, rotAngles=np.array([0, 0, 0]).all()):
    '''Duplicates unit cells to fill a given box'''
    print('hello I am working')

    def _rebuildAtom(atom, xyz):
        """Preserve original tuple shape and trailing meta data"""
        return tuple([atom[0], xyz[0], xyz[1], xyz[2], *list(atom[4:])])
    
    def _deduplicateAtomsPBC(atoms, latticeVec, tol=1e-3):
        '''Wrap all atoms to canonical unit cell and remove duplicates outside cell'''
        scale = 1.0/tol
        seen = {}
        unique = []

        for atom in atoms:
            key = (
                atom[0],
                int(round(atom[1] * scale)),
                int(round(atom[2] * scale)),
                int(round(atom[3] * scale)),
                tuple(atom[4:]),
            )
            if key not in seen:
                seen[key] = True
                unique.append(atom)

        return unique

    # Find the atom types in the tile
    atomNames = {atom[0] for atom in og}
    totalRadius = sum(radii[str.capitalize(name)] for name in atomNames)

    if not np.allclose(np.array(rotAngles), np.array([0, 0, 0])):
        latticeVec = calcLatticeVectors(cellDims[0], cellDims[1], cellDims[2], cellAngles[0], cellAngles[1], cellAngles[2])
        og, cellInfo = rotateUnitCell(latticeVec, og, rotAngles)
        cellDims = [cellInfo['a'], cellInfo['b'], cellInfo['c']]
        cellAngles = [cellInfo['alpha'], cellInfo['beta'], cellInfo['gamma']]

    latticeVec = calcLatticeVectors(cellDims[0], cellDims[1], cellDims[2], cellAngles[0], cellAngles[1], cellAngles[2])
    invLattice = np.linalg.inv(latticeVec)

    alpha = np.radians(cellAngles[0])
    beta = np.radians(cellAngles[1])
    gamma = np.radians(cellAngles[2])

    # Define the basis vectors
    a1 = latticeVec[0]
    a2 = latticeVec[1]
    a3 = latticeVec[2]

    # Wrap original atoms into on canonical unit cell
    wrappedOg = []
    wrapTol = 1e-8

    for atom in og:
        coord = np.array([atom[1], atom[2], atom[3]], dtype=float)
        
        # Cartesian to fraction for row based lattice vectors
        frac = coord @ invLattice

        # Wrap into [0,1]
        frac = frac % 1.0
        frac[np.isclose(frac, 1.0, atol=wrapTol)] = 0.0
        
        # Fractional to Cartesian
        wrappedCoord = frac @ latticeVec

        wrappedOg.append(_rebuildAtom(atom, wrappedCoord))
    canonicalOg = _deduplicateAtomsPBC(wrappedOg, latticeVec)

    # Calculate how many times to duplicate the unit cell in a given dimension
    dupeCount = [max(1, math.ceil(dims[i] / cellDims[i])) for i in range(3)]
    cellParams = f'CRYST1{dupeCount[0]*cellDims[0]:9.3f}{dupeCount[1]*cellDims[1]:9.3f}{dupeCount[2]*cellDims[2]:9.3f}  {cellAngles[0]:0.2f}  {cellAngles[1]:0.2f}  {cellAngles[2]:0.2f} P 1         1'

    filled = []

    for dx in range(dupeCount[0]):
        for dy in range(dupeCount[1]):
            for dz in range(dupeCount[2]):
                # Calculate the displacement for this duplication
                disp = dx*a1 + dy*a2 + dz*a3
                
                currentCell = []
                
                for atom in canonicalOg:
                    newCoord = np.array([atom[1], atom[2], atom[3]], dtype=float) + disp
                    shiftedCoord = np.array([
                        shape.xCoord + newCoord[0],
                        shape.yCoord + newCoord[1],
                        shape.zCoord + newCoord[2]
                    ], dtype=float)

                    # Check if the new atom fits within the box
                    currentCell.append(_rebuildAtom(atom, shiftedCoord))

                filled.append(currentCell)

    superLattice = np.array([
        dupeCount[0] * a1,
        dupeCount[1] * a2,
        dupeCount[2] * a3
    ])

    print("dims =", dims)
    print("cellDims =", cellDims)
    print("dupeCount =", dupeCount)
    print("a1 =", a1)
    print("a2 =", a2)
    print("a3 =", a3)
    print("superLattice =\n", superLattice)

    invSuper = np.linalg.inv(superLattice)

    wrappedFilled = []

    for cell in filled:
        wrappedCell = []

        for atom in cell:
            cart = np.array([
                atom[1] - shape.xCoord,
                atom[2] - shape.yCoord,
                atom[3] - shape.zCoord
            ], dtype=float)

            frac = cart @ invSuper
            frac = frac % 1.0
            frac[np.isclose(frac, 1.0, atol=wrapTol)] = 0.0

            wrappedCart = frac @ superLattice
            finalCoord = np.array([
                wrappedCart[0] + shape.xCoord,
                wrappedCart[1] + shape.yCoord,
                wrappedCart[2] + shape.zCoord
            ], dtype=float)

            wrappedCell.append(_rebuildAtom(atom, finalCoord))

        wrappedFilled.append(wrappedCell)
    return wrappedFilled, "molecule", cellParams

def unitCellSphere(shape, cellDims, cellAngles, og, radii, rotAngles=np.array([0, 0, 0]).all()):
    '''Duplicates unit cells to fill a given sphere'''

    # Box dimensions that completely contain the sphere
    boxDim = [2 * shape.radius] * 3

    if not np.allclose(np.array(rotAngles), np.array([0, 0, 0])):
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
    print(cellDims)
    print(cellAngles)
    boxDim = maxBound - minBound
    print(minBound, maxBound)
    if not np.allclose(np.array(rotAngles), np.array([0, 0, 0])):
        latticeVec = calcLatticeVectors(cellDims[0], cellDims[1], cellDims[2], cellAngles[0], cellAngles[1], cellAngles[2])
        og, cellInfo = rotateUnitCell(latticeVec, og, rotAngles)
        cellDims = [cellInfo['a'], cellInfo['b'], cellInfo['c']]
        print(cellDims)
        cellAngles = [cellInfo['alpha'], cellInfo['beta'], cellInfo['gamma']]
        print(cellAngles)

    alpha = np.radians(cellAngles[0])
    beta = np.radians(cellAngles[1])
    gamma = np.radians(cellAngles[2])

    dupeCount = [int(np.ceil(boxDim[i] / cellDims[i])) for i in range(3)]
    
    cellParams = f'CRYST1    {(dupeCount[0])*cellDims[0]: .3f}    {(dupeCount[1])*cellDims[1]: .3f}    {(dupeCount[2])*cellDims[2]: .3f}  {cellAngles[0]:0.2f}  {cellAngles[1]:0.2f}  {cellAngles[2]:0.2f} P1          1'

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
            for dz in range(dupeCount[2] + 5):
                disp = dx*a1 + dy*a2 + dz*a3 + minBound
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
