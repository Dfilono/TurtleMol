'''
This module handles a 3d box/cube to be filled
'''

class Box3d:
    '''
    Defines the 3D Box object
    '''
    def __init__(self, xCoord, yCoord, zCoord, dims):
        self.xCoord = xCoord
        self.yCoord = yCoord
        self.zCoord = zCoord
        self.length = dims[0]
        self.width = dims[1]
        self.height = dims[2]

    def volume(self):
        '''Returns the volume of a box'''
        return self.length * self.width * self.height

    def surfaceArea(self):
        '''Returns the surface area of a box'''
        face1 = 2 * self.length * self.width
        face2 = 2 * self.length * self.height
        face3 = 2 * self.width * self.height

        return face1 + face2 + face3

    def findCenter(self):
        '''Finds the center point of the Box'''
        corner1 = [self.xCoord, self.yCoord, self.zCoord]
        corner2 = [self.length, self.width, self.height]

        centerX = (corner1[0] + corner2[0]) / 2
        centerY = (corner1[1] + corner2[1]) / 2
        centerZ = (corner1[2] + corner2[2]) / 2

        return (centerX, centerY, centerZ)

    def origin(self):
        '''Find the origin'''
        return (self.xCoord, self.yCoord, self.zCoord)

def drawBox(params):
    '''
    Function sets the array that describes the dimensions of the box
    based on if the box is a cube or not
    '''
    if params['shape'].lower() == 'cube':
        return [float(params['sideLength']), float(params['sideLength']),
                float(params['sideLength'])]

    if params['shape'].lower() == 'box':
        return [float(params['Xlen']), float(params['Ylen']), float(params['Zlen'])]

    return 'ERROR: Shape not found\n'
