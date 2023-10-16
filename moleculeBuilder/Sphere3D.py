'''
This module handles a 3d sphere to be filled
'''

import math

class Sphere3d:
    '''
    Defines the 3D Sphere object
    '''
    def __init__(self, xCoord, yCoord, zCoord, radius):
        self.xCoord = xCoord
        self.yCoord = yCoord
        self.zCoord = zCoord
        self.radius = radius

    def volume(self):
        '''Returns the volume of a Sphere'''
        return (4/3) * math.pi * self.radius**3

    def surfaceArea(self):
        '''Returns the surface area of a Sphere'''
        return 4 * math.pi * self.radius**2

    def containsPoints(self, xCoord, yCoord, zCoord, atomRadius):
        '''Returns if the distance between points is greater than the radius'''
        distance = math.sqrt((xCoord - self.xCoord)**2 +
                             (yCoord - self.yCoord)**2 +
                             (zCoord - self.zCoord)**2)
        return distance <= (self.radius + atomRadius)

    def findCenter(self):
        '''Find the center of the sphere'''
        return (self.xCoord, self.yCoord, self.zCoord)

    def origin(self):
        '''Returns the origin of the sphere'''
        return (self.xCoord, self.yCoord, self.zCoord)
