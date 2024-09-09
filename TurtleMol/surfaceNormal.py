import numpy as np
from scipy.spatial.transform import Rotation as R
from .multMesh import computeCentroid

def alignToNormal(mesh, atoms):
    '''Aligns molecule's Z-Axis to the normal of the nearest face'''
    assert mesh is not None, "This function can only be used if the volume is a mesh object"

    centroid = computeCentroid(atoms)
    nearestPoint, normal = mesh.nearest.on_surface(centroid)
    zAxis = np.array([0, 0, 1])
    rotation, _ = R.align_vectors([normal], [zAxis])
    rotationMatrix = rotation.as_matrix()

    for atom in atoms:
        atomType = atom[0]
        points = atom[1:4]
        rotatedPoints = np.dot(points, rotationMatrix.T)
        rotatedPoints += centroid

        rotatedAtoms = [tuple(atomType) + tuple(rotatedPoints) if len(atom) == 4 else
                        tuple(atomType) + tuple(rotatedPoints) + tuple(atom[-1])]
        
    return rotatedAtoms

def placeOnSurfaceNormal(mesh, atoms):
    assert mesh is not None, "This function can only be used if the volume is a mesh object"

    centroid = computeCentroid(atoms)
    nearestPoint, normal = mesh.nearest.on_surface(centroid)
    zAxis = np.array([0, 0, 1])
    rotation, _ = R.align_vectors([normal], [zAxis])
    rotationMatrix = rotation.as_matrix()

    for i, atom in enumerate(atoms):
        atomType = atom[0]
        points = atom[1:4]
        rotatedPoints = np.dot(points - centroid, rotationMatrix.T)
        rotatedCentroid = np.mean(rotatedPoints, axis=0)

        translationToSurface = nearestPoint - rotatedCentroid

        newPos = rotatedPoints[i] + translationToSurface

        rotatedAtoms = [tuple(atomType) + tuple(newPos) if len(atom) == 4 else
                        tuple(atomType) + tuple(newPos) + tuple(atom[-1])]
        
    return rotatedAtoms

