'''Defines mesh'''

import numpy as np
from .readWriteFiles import readMesh

class mesh3D():
    '''Defines properties of a mesh'''

    def __init__(self, filePath):
        self.mesh = readMesh(filePath)

        # Mesh properties
        self.volume = self.mesh.volume
        self.surfaceArea = self.mesh.area
        self.isWaterTight = self.mesh.is_watertight

    def isInside(self, point):
        '''Check if point is inside the mesh'''
        return self.mesh.contains([point])[0]
