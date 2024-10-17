import numpy as np
import trimesh
from scipy.spatial.transform import Rotation as R
from .makeStruc import computeCentroid

def alignToNormal(mesh, atoms):
    '''Aligns molecule's Z-Axis to the normal of the nearest face'''
    assert mesh is not None, "This function can only be used if the volume is a mesh object"

    # Calculate the translation difference (how much the mesh has been moved)
    transDiff = mesh.bounds[0] - mesh.ogBounds[0]  # This is the amount the mesh was translated

    # Step 1: Temporarily translate the points back to the original mesh position
    revertedAtoms = []
    for atom in atoms:
        atomType = atom[0]
        points = np.array(atom[1:4])
        reverted_points = points - transDiff  # Reverse the translation
        revertedAtoms.append([atomType, reverted_points[0], reverted_points[1], reverted_points[2]] + list(atom[4:]))

    # Step 2: Compute the centroid of the translated points (in original position)
    centroid = np.array(computeCentroid(atoms))
    centroidShape = centroid.reshape((-1, 3))

    #normalIdx = trimesh.proximity.nearby_faces(mesh, centroidShape)
    nearestPoint, distance, triangleIdx = trimesh.proximity.closest_point(mesh, centroidShape)
    normals = mesh.face_normals
    zAxis = np.array([0, 0, 1])

    if np.allclose(zAxis, normals[triangleIdx]):
        rotationMatrix = np.eye(3)
    else:
        rotationMatrix = alignVectors(zAxis, normals[triangleIdx])
    
    rotated = []

    for atom in atoms:
        atomType = atom[0]
        points = atom[1:4]
        rotatedPoints = np.dot(points - centroid, rotationMatrix) + centroid

        if len(atom) == 4:
            rotated.append([atomType, rotatedPoints[0], rotatedPoints[1], rotatedPoints[2]])
        elif len(atom) == 5:
            rotated.append([atomType, rotatedPoints[0], rotatedPoints[1], rotatedPoints[2], atom[-1]])
    rotatedAtoms = [tuple(atom) for atom in rotated]
    return rotatedAtoms

def placeOnSurfaceNormal(mesh, atoms):
    '''Places the rotated molecule onto the nearest surface normal of the mesh'''

    # Calculate the translation difference (how much the mesh has been moved)
    transDiff = mesh.bounds[0] - mesh.ogBounds[0]  # This is the amount the mesh was translated

    # Step 1: Temporarily translate the points back to the original mesh position
    revertedAtoms = []
    for atom in atoms:
        atomType = atom[0]
        points = np.array(atom[1:4])
        reverted_points = points - transDiff  # Reverse the translation
        revertedAtoms.append([atomType, reverted_points[0], reverted_points[1], reverted_points[2]] + list(atom[4:]))

    # Step 2: Compute the centroid of the translated points (in original position)
    centroid = np.array(computeCentroid(revertedAtoms))
    centroidShape = centroid.reshape((-1, 3))

    # Find the nearest point and the normal of the nearest face on the original mesh
    nearestPoint, distance, triangleIdx = trimesh.proximity.closest_point(mesh, centroidShape)
    normals = mesh.face_normals
    zAxis = np.array([0, 0, 1])

    # Step 3: Rotate the molecule to align with the surface normal (still in original position)
    if np.allclose(zAxis, normals[triangleIdx]):
        rotationMatrix = np.eye(3)
    else:
        rotationMatrix = alignVectors(zAxis, normals[triangleIdx])

    rotatedAtoms = []
    for atom in revertedAtoms:
        atomType = atom[0]
        points = np.array(atom[1:4])
        rotatedPoints = np.dot(points - centroid, rotationMatrix) + centroid  # Rotate around the centroid
        rotatedAtoms.append([atomType, rotatedPoints[0], rotatedPoints[1], rotatedPoints[2]] + list(atom[4:]))

    # Step 4: Calculate the translation vector to place the rotated molecule on the surface (original position)
    rotatedCentroid = np.mean([atom[1:4] for atom in rotatedAtoms], axis=0)
    translationToSurface = nearestPoint - rotatedCentroid

    # Step 5: Apply the translation to the rotated atoms (still in original mesh position)
    translatedAtoms = []
    for atom in rotatedAtoms:
        atomType = atom[0]
        points = np.array(atom[1:4])
        newPos = points + translationToSurface
        translatedAtoms.append([atomType, newPos[0][0], newPos[0][1], newPos[0][2]] + list(atom[4:]))

    # Step 6: Translate both the nearest point and the molecule back to the current mesh position
    finalAtoms = []
    for atom in translatedAtoms:
        atomType = atom[0]
        points = np.array(atom[1:4])
        translatedPoints = points + transDiff  # Move the points back to the translated mesh position
        finalAtoms.append([atomType, translatedPoints[0], translatedPoints[1], translatedPoints[2]] + list(atom[4:]))

    # Return the final atoms placed on the translated mesh surface
    return [tuple(atom) for atom in finalAtoms]


def alignVectors(v1, v2):
    '''Generates a rotation matrix that aligns vector 1 to vector 2'''
    a, b = (v1 / np.linalg.norm(v1)).reshape(3), (v2 / np.linalg.norm(v2)).reshape(3)
    v = np.cross(a, b)
    c = np.dot(a, b)
    s = np.linalg.norm(v)

    # Handle the case when the vectors are parallel (s == 0)
    if np.isclose(s, 0):
        # Return identity matrix if no rotation is needed (vectors are parallel)
        return np.eye(3)
    
    kmat = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    rotationMatrix = np.eye(3) + kmat + kmat.dot(kmat) * ((1 - c) / (s**2))

    return rotationMatrix

