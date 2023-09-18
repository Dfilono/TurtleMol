import sys

'''
This module handles a 3d box/cube to be filled, and checks if atoms are still within the box
once they have been duplicated.
'''

class Box3d:
    def __init__(self, x, y, z, length, width, height):
        self.x = x
        self.y = y
        self.z = z
        self.length = length
        self.width = width
        self.height = height

    def is_inside(self, x, y, z):
        x_in_range = self.x <= x <= (self.x + self.length)
        y_in_range = self.y <= y <= (self.y + self.width)
        z_in_range = self.z <= z <= (self.z + self.height)

        if x_in_range and y_in_range and z_in_range:
            return True
        else:
            return False
        
def drawBox(params):
    if params['shape'].lower() == 'cube':
        return [float(params['sideLength']), float(params['sideLength']), float(params['sideLength'])]
    elif params['shape'].lower() == 'box':
        return [float(params['Xlen']), float(params['Ylen']), float(params['Zlen'])]
    else:
        print('Shape not found')
        sys.exit()