'''
These functions are used to repeat a pattern of a given molecule or molecules. 
The goal is to repeat the pattern exactly the specified number of times
within the constraints of a the defined shape. 
A 'Fill' option was also added so that one could just define a space,
and fill said space with provided structure as many times as it can fit.
'''

import math
from Box3D import Box3d
from Sphere3D import Sphere3d
from shiftBox import atomsFillBox, atomsRandBox, \
                     moleculesFillBox, moleculesRandBox
from shiftSphere import atomFillSphere, atomRandSphere, \
                        moleculeRandSphere, moleculeFillSphere
from setAtomProp import setAtomicRadius
from makeStruc import makeBase, calcNumMol

def drawMolBox(struc, tol, dims, maxattempts, numMol, baseStruc,
               randOrient, density, randFill):
    '''Utilized to place molecules in a box'''
    box = Box3d(0, 0, 0, dims)
    radii = setAtomicRadius()
    print(randFill)

    originalPoints = makeBase(struc)

    if density:
        numMol = calcNumMol(box, originalPoints, density)

    if not isinstance(numMol, int) and str(numMol).lower() != "fill":
        numMol = int(numMol)

    if randFill == 'False' or randFill == False:
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
                                          randOrient, numMol)
            return list(filledAtom), "molecule"

    elif randFill == 'True' or randFill is True:
        # Check if structure is monatomic or molecule
        if len(originalPoints) == 1:
            filledAtom = atomsRandBox(numMol, maxattempts,
                                          originalPoints, box, radii, tol)
            return list(filledAtom), "atom"

        if len(originalPoints) > 1:
            filledAtom = moleculesRandBox(numMol, maxattempts,
                                              originalPoints, box, radii, tol, baseStruc,
                                              randOrient)
            return list(filledAtom), "molecule"

    return "ERROR", "No atoms found"

def drawMolSphere(struc, tol, radius, center, maxattempts, numMol, baseStruc,
                  randOrient, density, randFill):
    '''Utilized to place molecules in a sphere'''
    sphere = Sphere3d(float(center[0]), float(center[1]), float(center[2]), radius)
    radii = setAtomicRadius()

    originalPoints = makeBase(struc)

    if density:
        numMol = calcNumMol(sphere, originalPoints, density)
        print(numMol)

    if not isinstance(numMol, int) and str(numMol).lower() != "fill":
        numMol = int(numMol)

    if randFill == 'False' or randFill is False:
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
                                            randOrient, numMol)
            return list(filledAtom), "molecule"

    elif randFill == 'True' or randFill == True:
        # Check if structure is monatomic or molecule
        if len(originalPoints) == 1:
            filledAtom = atomRandSphere(numMol, maxattempts,
                                            originalPoints, sphere, radii, tol)
            return list(filledAtom), "atom"

        if len(originalPoints) > 1:
            filledAtom = moleculeRandSphere(numMol, maxattempts,
                                                originalPoints, sphere, radii, tol, baseStruc,
                                                randOrient)
            return list(filledAtom), "molecule"

    return "ERROR: No atoms found"
