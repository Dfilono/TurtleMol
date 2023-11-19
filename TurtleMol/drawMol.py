'''
These functions are used to repeat a pattern of a given molecule or molecules. 
The goal is to repeat the pattern exactly the specified number of times
within the constraints of a the defined shape. 
A 'Fill' option was also added so that one could just define a space,
and fill said space with provided structure as many times as it can fit.
'''

import math
from .Box3D import Box3d, drawBox
from .Sphere3D import Sphere3d
from .shiftBox import atomsFillBox, atomsRandBox, \
                     moleculesFillBox, moleculesRandBox
from .shiftSphere import atomFillSphere, atomRandSphere, \
                        moleculeRandSphere, moleculeFillSphere
from .setAtomProp import setAtomicRadius
from .makeStruc import makeBase, calcNumMol, shiftPoints
from .shiftDensity import placeMols

def drawMolBox(struc, baseStruc, iparams):
    '''Utilized to place molecules in a box'''
    dims = drawBox(iparams)
    box = Box3d(0, 0, 0, dims)
    radii = setAtomicRadius(iparams['atomRadius'])
    numMol = iparams['numMolecules']
    tol = float(iparams['tol'])

    originalPoints = makeBase(struc)

    originalPoints = shiftPoints(originalPoints, box)

    if iparams['density']:
        filled, type = placeMols(box, originalPoints, iparams['density'],
                                 tol, "box", radii)
        return filled, type

    if not isinstance(numMol, int) and str(numMol).lower() != "fill":
        numMol = int(numMol)

    if iparams['randFill'] == 'False' or iparams['randFill'] is False:
        numXShifts = math.ceil(box.length / tol)
        numYShifts = math.ceil(box.height / tol)
        numZShifts = math.ceil(box.width / tol)

        # Check if structure is monatomic or molecule
        if len(originalPoints) == 1:
            filledAtom = atomsFillBox(numXShifts, numYShifts, numZShifts,
                                      tol, originalPoints, box, radii, numMol)
            return list(filledAtom), "atom"

        if len(originalPoints) > 1:
            filledAtom = moleculesFillBox(numXShifts, numYShifts, numZShifts,
                                          tol, originalPoints, box, radii, baseStruc,
                                          iparams['randomizeOrient'], numMol)
            return list(filledAtom), "molecule"

    elif iparams['randFill'] == 'True' or iparams['randFill'] is True:
        # Check if structure is monatomic or molecule
        if len(originalPoints) == 1:
            filledAtom = atomsRandBox(numMol, iparams['maxAttempts'],
                                          originalPoints, box, radii, tol)
            return list(filledAtom), "atom"

        if len(originalPoints) > 1:
            filledAtom = moleculesRandBox(numMol, iparams['maxAttempts'],
                                              originalPoints, box, radii, tol, baseStruc,
                                              iparams['randomizeOrient'])
            return list(filledAtom), "molecule"

    return "ERROR", "No atoms found"

def drawMolSphere(struc, baseStruc, iparams):
    '''Utilized to place molecules in a sphere'''
    sphere = Sphere3d(float(iparams['sphereCenter'][0]), float(iparams['sphereCenter'][1]),
                      float(iparams['sphereCenter'][2]), float(iparams['radius']))
    radii = setAtomicRadius(iparams['atomRadius'])
    numMol = iparams['numMolecules']
    tol = float(iparams['tol'])

    originalPoints = makeBase(struc)
    originalPoints = shiftPoints(originalPoints, sphere)

    if iparams['density']:
        filled, type = placeMols(sphere, originalPoints, iparams['density'],
                                 tol, "sphere", radii)
        return filled, type

    if not isinstance(numMol, int) and str(numMol).lower() != "fill":
        numMol = int(numMol)

    if iparams['randFill'] == 'False' or iparams['randFill'] is False:
        numShifts = math.ceil(2*sphere.radius / tol)

        # Check if structure is monatomic or molecule
        if len(originalPoints) == 1:
            filledAtom = atomFillSphere(numShifts, sphere,
                                        originalPoints, radii, tol,
                                        numMol)
            return list(filledAtom), "atom"

        if len(originalPoints) > 1:
            filledAtom = moleculeFillSphere(numShifts, sphere,
                                            originalPoints, radii, tol, baseStruc,
                                            iparams['randomizeOrient'], numMol)
            return list(filledAtom), "molecule"

    elif iparams['randFill'] == 'True' or iparams['randFill'] is True:
        # Check if structure is monatomic or molecule
        if len(originalPoints) == 1:
            filledAtom = atomRandSphere(numMol, iparams['maxAttempts'],
                                            originalPoints, sphere, radii, tol)
            return list(filledAtom), "atom"

        if len(originalPoints) > 1:
            filledAtom = moleculeRandSphere(numMol, iparams['maxAttempts'],
                                                originalPoints, sphere, radii, tol, baseStruc,
                                                iparams['randomizeOrient'])
            return list(filledAtom), "molecule"

    return "ERROR: No atoms found"
