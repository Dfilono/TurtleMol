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
    #print(normals[normalIdx[0][0]])
    zAxis = np.array([0, 0.5, 0.5])
    rotationMatrix = alignVectors([zAxis], [normals[normalIdx[0][0]]])
    
    
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
    #print(len(rotated))
    rotatedAtoms = [tuple(atom) for atom in rotated]
    #print(len(rotatedAtoms))
    return rotatedAtoms

def placeOnSurfaceNormal(mesh, atoms):
    assert mesh is not None, "This function can only be used if the volume is a mesh object"

    centroid = np.array(computeCentroid(atoms))
    centroidShape = centroid.reshape((-1, 3))
    normalIdx = trimesh.proximity.nearby_faces(mesh, centroidShape)
    nearestPoint, distance, triangleIdx = trimesh.proximity.closest_point(mesh, centroidShape)
    normals = mesh.face_normals
    zAxis = np.array([0, 0.5, 0.5])
    rotationMatrix = alignVectors([zAxis], [normals[normalIdx[0][0]]])

    rotatedAtoms = []

    for i, atom in enumerate(atoms):
        atomType = atom[0]
        points = np.array(atom[1:4])
        rotatedPoints = rotationMatrix.dot(points - centroid)
        if len(atom) == 4:
            rotatedAtoms.append([atomType, *rotatedPoints])
        elif len(atom) == 5:
            rotatedAtoms.append([atomType, *rotatedPoints, atom[-1]])

    rotatedCentroid = np.mean([atom[1:4] for atom in rotatedAtoms], axis=0)

    translationToSurface = nearestPoint - rotatedCentroid

    translated = []

    for atom in rotatedAtoms:
        atomType = atom[0]
        points = np.array(atom[1:4])
        newPos = points + translationToSurface
        if len(atom) == 4:
            translated.append([atomType, newPos[0][0], newPos[0][1], newPos[0][2]])
        elif len(atom) == 5:
            translated.append([atomType, newPos[0][0], newPos[0][1], newPos[0][2], atom[-1]])

    translatedAtoms = [tuple(atom) for atom in translated]
        
    return translatedAtoms

def alignVectors(v1, v2):
    '''Generates a rotation matrix that aligns vector 1 to vector 2'''
    a, b = (v1 / np.linalg.norm(v1)).reshape(3), (v2 / np.linalg.norm(v2)).reshape(3)
    v = np.cross(a, b)
    c = np.dot(a, b)
    s = np.linalg.norm(v)
    #print(s)
    kmat = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    rotationMatrix = np.eye(3) + kmat + kmat.dot(kmat) * ((1 - c) / (s**2))

    return rotationMatrix

