'''
This module handles a 3d box/cube to be filled, and checks if atoms are still within the box
once they have been duplicated.
'''

import sys

class Box3d:
    def __init__(self, x, y, z, length, width, height):
        self.x = x
        self.y = y
        self.z = z
        self.length = length
        self.width = width
        self.height = height

'''
Function sets the array that describes the dimensions of the box based on if the box is a cube or not
'''
def drawBox(params):
    if params['shape'].lower() == 'cube':
        return [float(params['sideLength']), float(params['sideLength']), 
                float(params['sideLength'])]

    if params['shape'].lower() == 'box':
        return [float(params['Xlen']), float(params['Ylen']), float(params['Zlen'])]