'''Sets up the structure files'''
import random
import math
import numpy as np
from .setAtomProp import setAtomicMass

def makeBase(baseStruc):
    '''Convert dataframe to list of points'''
    atomBase = baseStruc['Atom'].values.tolist()
    xBase = baseStruc['X'].values.tolist()
    yBase = baseStruc['Y'].values.tolist()
    zBase = baseStruc['Z'].values.tolist()

    if 'Residue' not in baseStruc.columns:
        return [(atomBase[i], xBase[i], yBase[i], zBase[i]) for i in range(len(atomBase))]

    residue = baseStruc['Residue'].values.tolist()
    return [(atomBase[i], xBase[i], yBase[i], zBase[i], residue[i]) for i in range(len(atomBase))]

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
    newCoords = [(coord[0], coord[1] + displacement[0], coord[2] + displacement[1],
                coord[3] + displacement[2]) for coord in struc]

    return newCoords

def shiftPoints(points, shape):
    '''Sifts points to the bottom-left front corner of a 3D shape'''
    pointNames = [point[0] for point in points]
    coords = np.array([point[1:4] for point in points], dtype=float)

    if len(points[0]) == 5:
        pointRes = [point[4] for point in points]
    else:
        pointRes = None

    # Determine the bounds of the points
    minXPoints = np.min(coords[:, 0])
    minYPoints = np.min(coords[:, 1])
    minZPoints = np.min(coords[:, 2])

    # Calculate the translation required
    translation = np.array([
        shape.origin()[0] - minXPoints,
        shape.origin()[1] - minYPoints,
        shape.origin()[2] - minZPoints,
    ])

    # Shift the points
    shiftedCoords = coords + translation
    if pointRes:
        shiftedPoints = [[name] + coord.tolist() + [residue]
                         for name, coord, residue in
                         zip(pointNames, shiftedCoords, pointRes)]
    else:
        shiftedPoints = [[name] + coord.tolist()
                         for name, coord in zip(pointNames, shiftedCoords)]

    return shiftedPoints

def randReorient(mol): # NOTE Currently does not work and is disabled
    '''Randomly reorients a molecule around geometric center'''
    if len(mol) == 0:
        return mol

    atomNames = [atom[0] for atom in mol]
    atomInfo = [atom[4:] if len(atom) > 4 else [] for atom in mol]
    points = np.array([atom[1:4] for atom in mol], dtype=np.float64)
    centroid = points.mean(axis=0) # Calculate center of points
    points -= centroid # Translate the points to origin

    # Perform the rotaton
    # Generate random rotation angles in radians
    thetaX = random.uniform(0, 2*np.pi)
    thetaY = random.uniform(0, 2*np.pi)
    thetaZ = random.uniform(0, 2*np.pi)

    # Create the rotation matrix for the x-axis
    rx = np.array([[1, 0, 0],
                  [0, np.cos(thetaX), -np.sin(thetaX)],
                  [0, np.sin(thetaX), np.cos(thetaX)]],
                  dtype=np.float64)

    # Create rotation matrix for the y-axis
    ry = np.array([[np.cos(thetaY), 0, np.sin(thetaY)],
                  [0, 1, 0],
                  [-np.sin(thetaY), 0, np.cos(thetaY)]],
                  dtype=np.float64)

    # Create rotation matrix for the z-axis
    rz = np.array([[np.cos(thetaZ), -np.sin(thetaZ), 0],
                  [np.sin(thetaZ), np.cos(thetaZ), 0],
                  [0, 0, 1]],
                  dtype=np.float64)

    # Combined rotation matrix
    R = rz @ ry @ rx

    # Apply the rotation matrix to the points
    rotatedPoints = np.dot(points, R.T)

    # Translate the points back to position
    rotatedPoints += centroid

    rotatedAtoms = [tuple(name) + tuple(coord.tolist()) + tuple(info)
                    for name, coord, info in zip(atomNames, rotatedPoints, atomInfo)]

    rotatedAtoms = [tuple(atom) for atom in rotatedAtoms]

    return rotatedAtoms

def calcDensity(shape, mol):
    '''Calculates the density of a given structure'''
    vol = shape.volume() * 1e-24 # mL
    mass = 0
    atomicMass = setAtomicMass()  # g

    for atom in mol:
        mass += atomicMass[atom[0]]
    
    moles = mass # g

    return mass/vol # g/mL

def calcNumMol(shape, mol, density):
    '''
    Calulates the number of molecules needed in a box to match
    the defined density based off of the defined volume
    '''
    vol = shape.volume() # A^3
    rho = density # g/mL
    rho *= 1e-24 # g/A^3
    atomicMass = setAtomicMass() # g/mol
    molarMass = 0 
    for atom in mol:
        molarMass += atomicMass[atom[0]]
    
    numMol = (rho * 6.022e23 * vol)/molarMass # molecules

    return int(numMol)

def calcDistance(shape, mol, density, shapeType):
    '''
    Calculates the distance between molecules for the correct density
    '''
    if str(shapeType).lower() == "box":
        x, y, z = shape.length, shape.width, shape.height

        numMol = calcNumMol(shape, mol, density)

        # Calculate the molecules to place in each dimension
        pointsPerDimension = round(numMol ** (1/3))

        # Generate grid points
        XCount = np.linspace(0, x, pointsPerDimension, endpoint=False)
        YCount = np.linspace(0, y, pointsPerDimension, endpoint=False)
        ZCount = np.linspace(0, z, pointsPerDimension, endpoint=False)

        # Create a meshgrid and reshape it to get coordiantes
        xx, yy, zz = np.meshgrid(XCount, YCount, ZCount, indexing='ij')
        coords = np.vstack([xx.ravel(), yy.ravel(), zz.ravel()]).T

        return coords[:numMol].tolist()
    
    elif str(shapeType).lower() == "sphere":

        def insideSphere(x, y, z, radius):
            return x**2 + y**2 + z**2 <= radius**2

        r = shape.radius

        numMol = calcNumMol(shape, mol, density)

        # Create cubic geometry to enclose the sphere
        cubeSide = 2*r

        # Calculate total points in the cube
        rho = numMol / shape.volume()
        totalPoints = int(rho * cubeSide**3)

        # Generate points within cube
        XCount = np.linspace(-r, r, int(totalPoints**(1/3)))
        YCount = np.linspace(-r, r, int(totalPoints**(1/3)))
        ZCount = np.linspace(-r, r, int(totalPoints**(1/3)))

        # Create meshgrid and reshape to get coordinates
        xx, yy, zz = np.meshgrid(XCount, YCount, ZCount, indexing='ij')
        coords = np.vstack([xx.ravel(), yy.ravel(), zz.ravel()]).T

        # Find points inside sphere
        insidePoints = [point for point in coords if insideSphere(*point, r)]
        
        return insidePoints[:numMol].tolist()