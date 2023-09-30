'''
This module handles a 3d box/cube to be filled, and checks if atoms are still within the box
once they have been duplicated.
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

def drawBox(params):
    '''
    Function sets the array that describes the dimensions of the box based on if the box is a cube or not
    '''
    if params['shape'].lower() == 'cube':
        return [float(params['sideLength']), float(params['sideLength']),
                float(params['sideLength'])]

    if params['shape'].lower() == 'box':
        return [float(params['Xlen']), float(params['Ylen']), float(params['Zlen'])]
