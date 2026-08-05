'''This module is for replicating a unit cell if it is defined in a given pdb'''

import numpy as np
from scipy.spatial import cKDTree
from collections import defaultdict
import math
import trimesh
from .shiftBox import inBox
from .treeManager import PendingKDManager
from .makeStruc import calcLatticeVectors, rotateUnitCell

# Box

def unitCellBox(shape, dims, cellDims, cellAngles, og, radii, rotAngles=np.array([0, 0, 0]).all(), conect=None):
    '''Duplicates unit cells to fill a given box'''

    def _rebuildAtom(atom, xyz):
        """Preserve original tuple shape and trailing meta data"""
        return tuple([atom[0], xyz[0], xyz[1], xyz[2], *list(atom[4:])])
    
    def _deduplicateAtomsPBC(atoms, latticeVec, tol=1e-3):
        """
        Deduplicate atoms under periodic boundary conditions using
        minimum-image distance in fractional coordinates.
        """
        coords = np.array([[a[1], a[2], a[3]] for a in atoms], dtype=float)
        invLattice = np.linalg.inv(latticeVec)

        # Fractional coords in [0,1)
        frac = coords @ invLattice
        frac = frac % 1.0

        n = len(atoms)
        parent = list(range(n))

        def find(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x

        def union(a, b):
            ra, rb = find(a), find(b)
            if ra != rb:
                parent[rb] = ra

        # Brute force is safest for correctness here.
        # If needed, this can be optimized later.
        for i in range(n):
            for j in range(i + 1, n):
                if atoms[i][0] != atoms[j][0]:
                    continue

                df = frac[i] - frac[j]
                df -= np.round(df)          # minimum-image fractional delta
                dr = df @ latticeVec        # back to Cartesian
                dist = np.linalg.norm(dr)

                if dist < tol:
                    union(i, j)

        unique = []
        kept_indices = []
        seen = set()

        for i in range(n):
            r = find(i)
            if r not in seen:
                seen.add(r)
                unique.append(atoms[i])
                kept_indices.append(i)

        return unique, kept_indices
    
    def _filter_intra_image_conect(conect, kept_indices):
        """
        Keep only bonds where both atoms survived canonical wrapping/dedup.
        Assumes conect keys are 1-based indices matching original og order.
        Returns:
            baseConect using 1-based indices into canonicalOg
        """
        if not conect:
            return {}

        old_to_new = {old_idx + 1: new_idx + 1 for new_idx, old_idx in enumerate(kept_indices)}
        baseConect = defaultdict(set)

        for src_old, nbrs in conect.items():
            if src_old not in old_to_new:
                continue
            src_new = old_to_new[src_old]

            for dst_old in nbrs:
                if dst_old not in old_to_new:
                    continue
                dst_new = old_to_new[dst_old]

                if src_new == dst_new:
                    continue

                baseConect[src_new].add(dst_new)
                baseConect[dst_new].add(src_new)

        return dict(baseConect)

    def _replicate_conect(baseConect, atoms_per_cell, n_cells):
        """
        Replicate intra-image connectivity into each duplicated cell.
        baseConect uses 1-based indices inside one canonical cell.
        """
        if not baseConect:
            return {}

        newConect = defaultdict(set)

        for cell_idx in range(n_cells):
            offset = cell_idx * atoms_per_cell
            for src, nbrs in baseConect.items():
                src_new = src + offset
                for dst in nbrs:
                    dst_new = dst + offset
                    newConect[src_new].add(dst_new)
                    newConect[dst_new].add(src_new)

        return dict(newConect)

    # Find the atom types in the tile
    atomNames = {atom[0] for atom in og}
    totalRadius = sum(radii[str.capitalize(name)] for name in atomNames)

    if not np.allclose(np.array(rotAngles), np.array([0, 0, 0])):
        latticeVec = calcLatticeVectors(
            cellDims[0], cellDims[1], cellDims[2],
            cellAngles[0], cellAngles[1], cellAngles[2]
        )
        og, cellInfo = rotateUnitCell(latticeVec, og, rotAngles)
        cellDims = [cellInfo['a'], cellInfo['b'], cellInfo['c']]
        cellAngles = [cellInfo['alpha'], cellInfo['beta'], cellInfo['gamma']]

    latticeVec = calcLatticeVectors(
        cellDims[0], cellDims[1], cellDims[2],
        cellAngles[0], cellAngles[1], cellAngles[2]
    )
    invLattice = np.linalg.inv(latticeVec)

    a1 = latticeVec[0]
    a2 = latticeVec[1]
    a3 = latticeVec[2]

    # Wrap original atoms into one canonical unit cell
    wrappedOg = []
    wrapTol = 1e-8

    for atom in og:
        coord = np.array([atom[1], atom[2], atom[3]], dtype=float)

        # Cartesian -> fractional for row-based lattice vectors
        frac = coord @ invLattice

        # Wrap into [0,1)
        frac = frac % 1.0
        frac[np.isclose(frac, 1.0, atol=wrapTol)] = 0.0

        # Fractional -> Cartesian
        wrappedCoord = frac @ latticeVec

        wrappedOg.append(_rebuildAtom(atom, wrappedCoord))

    # Deduplicate wrapped canonical cell
    canonicalOg, kept_indices = _deduplicateAtomsPBC(wrappedOg, latticeVec, tol=0.15)

    # Preserve only intra-image conects among surviving canonical atoms
    baseConect = _filter_intra_image_conect(conect, kept_indices)

    # Use ceil so the requested box is actually covered
    dupeCount = [max(1, math.ceil(dims[i] / cellDims[i])) for i in range(3)]

    cellParams = (
        f'CRYST1{dupeCount[0]*cellDims[0]:9.3f}'
        f'{dupeCount[1]*cellDims[1]:9.3f}'
        f'{dupeCount[2]*cellDims[2]:9.3f}  '
        f'{cellAngles[0]:0.2f}  {cellAngles[1]:0.2f}  {cellAngles[2]:0.2f} P 1         1'
    )

    filled = []

    for dx in range(dupeCount[0]):
        for dy in range(dupeCount[1]):
            for dz in range(dupeCount[2]):
                disp = dx * a1 + dy * a2 + dz * a3
                currentCell = []

                for atom in canonicalOg:
                    newCoord = np.array([atom[1], atom[2], atom[3]], dtype=float) + disp
                    shiftedCoord = np.array([
                        shape.xCoord + newCoord[0],
                        shape.yCoord + newCoord[1],
                        shape.zCoord + newCoord[2]
                    ], dtype=float)

                    currentCell.append(_rebuildAtom(atom, shiftedCoord))

                filled.append(currentCell)

    # Wrap all atoms into final supercell box
    superLattice = np.array([
        dupeCount[0] * a1,
        dupeCount[1] * a2,
        dupeCount[2] * a3
    ])
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


    atoms_per_cell = len(canonicalOg)

    newConect = _replicate_conect(baseConect, atoms_per_cell, len(wrappedFilled))

    if conect is None:
        return wrappedFilled, "molecule", cellParams
    return wrappedFilled, "molecule", cellParams, newConect

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
