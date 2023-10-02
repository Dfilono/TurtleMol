'''Sets up the structure files'''

def makeBase(baseStruc):
    '''Convert dataframe to list of points'''
    atomBase = baseStruc['Atom'].values.tolist()
    xBase = baseStruc['X'].values.tolist()
    yBase = baseStruc['Y'].values.tolist()
    zBase = baseStruc['Z'].values.tolist()

    return [(atomBase[i], xBase[i], yBase[i], zBase[i]) for i in range(len(atomBase))]

def calcCenter(coords):
    '''Calculate geometric center'''
    if not coords:
        return None

    numCoords = len(coords)
    centerX = sum(coord[1] for coord in coords) / numCoords
    centerY = sum(coord[2] for coord in coords) / numCoords
    centerZ = sum(coord[3] for coord in coords) / numCoords

    return (centerX, centerY, centerZ)    

def reCenter(struc, shape):
    '''Set center structure coordiantes to center of shape'''
    currentCenter = calcCenter(struc)
    shapeCenter = shape.findCenter()
    displacement = (
        shapeCenter[0] - currentCenter[0],
        shapeCenter[1] - currentCenter[1],
        shapeCenter[2] - currentCenter[2]
    )
    newCoords = [(coord[0], coord[1] + displacement[0], coord[2] + displacement[1], coord[3] + displacement[2]) for coord in struc]

    return newCoords