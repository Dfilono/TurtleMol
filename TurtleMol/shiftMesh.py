'''Fills an arbitrary mesh with molecules'''

import random
import numpy as np
from .makeStruc import makeBase, reCenter, Reorient
from .surfaceNormal import placeOnSurfaceNormal, alignToNormal
from. treeManager import PendingKDManager

def atomsFillMesh(mesh, og, tol, radii, numMol):
    '''Fills mesh with single atoms'''
    filled = []

    if str(numMol).lower() == 'fill':
        numMol = 10**100

    # Create KD-tree for filledAtoms
    manager = PendingKDManager(radii, rebuildRate=500)

    # Determine bounds of mesh
    bounds = mesh.bounds
    minBound, maxBound = bounds[0], bounds[1]
    print(minBound, maxBound)

    # Generate grid of points
    xs = np.arange(minBound[0], maxBound[0], tol)
    ys = np.arange(minBound[1], maxBound[1], tol)
    zs = np.arange(minBound[2], maxBound[2], tol)

    X, Y, Z = np.meshgrid(xs, ys, zs, indexing="ij")
    points = np.stack([X.ravel(), Y.ravel(), Z.ravel()], axis=1).astype(np.float64)

    # Vectorized inside test
    # Depending on mesh wrapper
    # if mesh.meshBound is a trimesh object: mesh.meshBound.contains(points)
    # if mesh itself is trimesh: mesh.contains(points)
    mask = mesh.isInsideMany(points)
    insidePoints = points[mask]
    if insidePoints.shape[0] == 0:
        return filled
    
    atomTemplate = og[0] # pick the first atom type as the one we're placing
    atomType = atomTemplate[0]
    extra = atomTemplate[4:] # handles optional fields

    # Check each point in the grid
    for (x, y, z) in insidePoints:
        if len(filled) >= numMol:
            break

        newAtom = (atomType, float(x), float(y), float(z), *extra)
        newMol = [newAtom]

        if not manager.overlapsMolecule(newMol, tol):
            filled.append(newMol)
            manager.addMolecule(newMol)
            manager.maybeRebuild()
    
    return filled

def moleculesFillMesh(mesh, og, tol, radii, numMol, baseStruc, 
                      randOrient, rotAngles, alignNormal=False, onSurface=False):
    '''Fills mesh with molecules'''
    filled = []

    if str(numMol).lower() == 'fill':
        numMol = 10**100

    # Use pending-buffer KD manager
    manager = PendingKDManager(radii, rebuildRate=500)

    if baseStruc is not None:
        base = makeBase(baseStruc)
        filled.append(reCenter(base, mesh))
        manager.addMolecule(base)
        manager.maybeRebuild(force=True) #ensure the base is comitted

    # Determine bounds of mesh
    bounds = mesh.bounds
    minBound, maxBound = bounds[0], bounds[1]

    xs = np.arange(minBound[0], maxBound[0], tol)
    ys = np.arange(minBound[1], maxBound[1], tol)
    zs = np.arange(minBound[2], maxBound[2], tol)

    X, Y, Z = np.meshgrid(xs, ys, zs, indexing="ij")
    anchors = np.stack([X.ravel(), Y.ravel(), Z.ravel()], axis=1).astype(np.float64)

    ogXYZ = np.array([[a[1], a[2], a[3]] for a in og], dtype=float)
    center = ogXYZ.mean(axis=0)

    ogRel = []
    for a in og:
        if len(a) == 4:
            ogRel.append((a[0], float(a[1] - center[0]), float(a[2] - center[1]), float(a[3] - center[2])))
        else:
            ogRel.append((a[0], float(a[1]-center[0]), float(a[2]-center[1]), float(a[3]-center[2]), a[4]))
    
    offsets = np.array([[a[1], a[2], a[3]] for a in ogRel], dtype=np.float64)
    atomTypes = [a[0] for a in ogRel]
    extras = [a[4] if len(a) == 5 else None for a in ogRel]

    for ax, ay, az in anchors:
        if len(filled) >= numMol:
            break

        points = offsets + np.array([ax, ay, az], dtype=np.float64)

        if not np.all(mesh.isInsideMany(points)):
            continue

        newMol = []
        for (x, y, z), t, ex in zip(points, atomTypes, extras):
            if ex is None:
                newMol.append((t, float(x), float(y), float(z)))
            else:
                newMol.append((t, float(x), float(y), float(z), ex))
        if randOrient:
            newMol = Reorient(newMol, randRotate=True)
        if not np.allclose(np.array(rotAngles), np.array([0,0,0])):
            newMol = Reorient(newMol, angles=rotAngles)
        if alignNormal:
            newMol = alignToNormal(mesh, newMol)
        if onSurface:
            newMol = placeOnSurfaceNormal(mesh, newMol)

        if not manager.overlapsMolecule(newMol, tol):
            filled.append(newMol)
            manager.addMolecule(newMol)
            manager.maybeRebuild()
    
    return filled

def atomsRandMesh(mesh, og, tol, radii, numMol, maxAttempts):
    '''Randomly places atoms in a mesh'''
    filled = []
    attempts = 0

    # Use pending-buffer KD manager
    manager = PendingKDManager(radii, rebuildRate=500)

    # Determine bounds of mesh
    bounds = mesh.bounds
    minBound, maxBound = bounds[0], bounds[1]

    while len(filled) < numMol and attempts <= maxAttempts:
        newMol = []
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

                newMol.append(atomData)
                if (numMol > len(filled) and not manager.overlapsMolecule(newMol, tol)):
                    filled.append(newMol)

                    # Rebuild KDTree with newly added atoms
                    manager.addMolecule(newMol)
                    manager.maybeRebuild() 

                if len(filled) >= numMol:
                    return filled

        attempts += 1
    return filled

def moleculesRandMesh(mesh, og, tol, radii, numMol, baseStruc,
                      randOrient, rotAngles, maxAttempts):
    '''Randomly places molecules in a mesh'''
    filled = []
    attempts = 0

    # Use pending-buffer KD manager
    manager = PendingKDManager(radii, rebuildRate=500)

    if baseStruc is not None:
        base = makeBase(baseStruc)
        filled.append(reCenter(base, mesh))
        manager.addMolecule(base)
        manager.maybeRebuild(force=True) #ensure the base is comitted

    # ensure og is relative offsets around 0
    ogXYZ = np.array([[a[1], a[2], a[3]] for a in og], dtype=np.float64)
    center = ogXYZ.mean(axis=0)

    ogRel = []
    for a in og:
        if len(a) == 4:
            ogRel.append((a[0], float(a[1]-center[0]), float(a[2]-center[1]), float(a[3]-center[2])))
        else:
            ogRel.append((a[0], float(a[1]-center[0]), float(a[2]-center[1]), float(a[3]-center[2]), a[4]))

    # Determine bounds of mesh
    bounds = mesh.bounds
    minBound, maxBound = bounds[0], bounds[1]

    offsets = np.array([[a[1], a[2], a[3]] for a in ogRel], dtype=np.float64)
    atomTypes = [a[0] for a in ogRel]
    extras = [a[4] if len(a) == 5 else None for a in ogRel]

    while len(filled) < numMol and attempts <= maxAttempts:
        newMol = []

        ax = random.uniform(minBound[0], maxBound[0])
        ay = random.uniform(minBound[1], maxBound[1])
        az = random.uniform(minBound[2], maxBound[2])

        points = offsets + np.array([ax, ay, az], dtype=np.float64)

        if not np.all(mesh.isInsideMany(points)):
            attempts += 1
            continue

        for (x, y, z), t, ex in zip(points, atomTypes, extras):
            if ex is None:
                newMol.append((t, float(x), float(y), float(z)))
            else:
                newMol.append((t, float(x), float(y), float(z), ex))

        if randOrient and len(newMol) == len(og):
            newMol = Reorient(newMol, randRotate=True)

        if not np.allclose(np.array(rotAngles), np.array([0, 0, 0])) and len(newMol) == len(og):
            newMol = Reorient(newMol, angles=rotAngles)

        if (numMol > len(filled) and not manager.overlapsMolecule(newMol, tol)):
            filled.append(newMol)

            # Rebuild KDTree with newly added atoms
            manager.addMolecule(newMol)
            manager.maybeRebuild()

        if len(filled) >= numMol:
            return filled

        attempts += 1
    return list(filled)
