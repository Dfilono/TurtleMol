import numpy as np
import trimesh
from scipy.spatial.transform import Rotation as R
from .makeStruc import computeCentroid

def alignToNormal(mesh, atoms):
    '''Aligns molecule's Z-Axis to the normal of the nearest face'''
    assert mesh is not None, "This function can only be used if the volume is a mesh object"

    centroid = np.array(computeCentroid(atoms))
    centroidShape = centroid.reshape((-1, 3))
    normalIdx = trimesh.proximity.nearby_faces(mesh, centroidShape)
    normals = mesh.face_normals
    zAxis = np.array([0, 0, 1])
    rotation, _ = R.align_vectors([normals[normalIdx[0][0]]], [zAxis])
    rotationMatrix = rotation.as_matrix()
    
    rotated = []

    for atom in atoms:
        atomType = atom[0]
        points = atom[1:4]
        rotatedPoints = np.dot(points, rotationMatrix)
        rotatedPoints += centroid

        if len(atom) == 4:
            rotated.append([atomType, rotatedPoints[0], rotatedPoints[1], rotatedPoints[2]])
        elif len(atom) == 5:
            rotated.append([atomType, rotatedPoints[0], rotatedPoints[1], rotatedPoints[2], atom[-1]])

    rotatedAtoms = [tuple(atom) for atom in rotated]
        
    return rotatedAtoms

def placeOnSurfaceNormal(mesh, atoms):
    assert mesh is not None, "This function can only be used if the volume is a mesh object"

    centroid = np.array(computeCentroid(atoms))
    centroidShape = centroid.reshape((-1, 3))
    normalIdx = trimesh.proximity.nearby_faces(mesh, centroidShape)
    nearestPoint, distance, triangleIdx = trimesh.proximity.closest_point(mesh, centroidShape)
    normals = mesh.face_normals
    zAxis = np.array([0, 0, 1])
    rotation, _ = R.align_vectors([normals[normalIdx[0][0]]], [zAxis])
    rotationMatrix = rotation.as_matrix()

    rotated = []

    for i, atom in enumerate(atoms):
        atomType = atom[0]
        points = atom[1:4]

        rotatedPoints = rotationMatrix.dot(points - centroid)
        #rotatedPoints = np.dot(points - centroid, rotationMatrix)
        rotatedCentroid = np.mean(rotatedPoints, axis=0)

        translationToSurface = nearestPoint - rotatedCentroid

        newPos = rotatedPoints[i] + translationToSurface

        if len(atom) == 4:
            rotated.append([atomType, newPos[0][0], newPos[0][1], newPos[0][2]])
        elif len(atom) == 5:
            rotated.append([atomType, newPos[0][0], newPos[0][1], newPos[0][2], atom[-1]])

    rotatedAtoms = [tuple(atom) for atom in rotated]
        
    return rotatedAtoms

