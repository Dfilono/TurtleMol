'''Defines mesh'''

import numpy as np
from .readWriteFiles import readMesh
import trimesh

class mesh3D():
    '''Defines properties of a mesh'''

    def __init__(self, filePath, scale):
        self.mesh = readMesh(filePath)
        self.scale = float(scale)
        self.mesh = self.mesh.apply_scale(self.scale)

        # Mesh properties
        self.isWaterTight = self.mesh.is_watertight
        self.bounds = self.mesh.bounds
        self.meshOrigin = self.origin()

        self.translate()

        # Scaled mesh
        self.meshBound = self.mesh.apply_scale((1.1, 1.1, 1))

    def isInside(self, point):
        '''Check if point is inside the mesh or on the boundary'''
        # Checks if the point is inside the mesh
        inside = self.meshBound.contains([point])[0]

        # Check if the point is on the boundary of the mesh
        # We will get the closest point on the mesh to our point and check if they are the same
        closestPoint, distance, _ = trimesh.proximity.closest_point(self.mesh, [point])
        onBound = np.isclose(distance[0], 0, atol=1e-3)

        return inside or onBound
    
    def volume(self):
        return self.mesh.volume
    
    def surfaceArea(self):
        return self.mesh.area
    
    def origin(self):
        return [0, 0, 0]
    
    def findCenter(self):
        return self.origin
    
    def translate(self):
        self.mesh.apply_translation(-self.bounds[0])
        self.bounds = self.mesh.bounds