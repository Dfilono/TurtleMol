'''
These functions are used to repeat a pattern of a given molecule or molecules. 
The goal is to repeat the pattern exactly the specified number of times
within the constraints of a the defined shape. 
A 'Fill' option was also added so that one could just define a space,
and fill said space with provided structure as many times as it can fit.
'''

import math
import pandas as pd
from Box3D import Box3d
from Sphere3D import Sphere3d
from shiftBox import atomsFillBox, atomsDefiniteBox, \
                     moleculesFillBox, moleculesDefiniteBox
from shiftSphere import atomFillSphere, atomDefiniteSphere, \
                        moleculeDefiniteSphere, moleculeFillSphere
from setAtomProp import setAtomicRadius
from makeStruc import makeBase

def drawMolBox(struc, tol, dims, maxattempts, numMol, baseStruc, randOrient):
    '''Utilized to place molecules in a box'''
    box = Box3d(0, 0, 0, dims)
    radii = setAtomicRadius()

    originalPoints = makeBase(struc)

    if not isinstance(numMol, int) and str(numMol).lower() != "fill":
        numMol = int(numMol)

    if str(numMol).lower() == "fill":
        numXShifts = math.ceil(box.length / tol)
        numYShifts = math.ceil(box.height / tol)
        numZShifts = math.ceil(box.width / tol)

        # Check if structure is monatomic or molecule
        if len(originalPoints) == 1:
            filledAtom = atomsFillBox(numXShifts, numYShifts, numZShifts,
                                      tol, originalPoints, box, radii)
            return list(filledAtom), "atom"

        if len(originalPoints) > 1:
            filledAtom = moleculesFillBox(numXShifts, numYShifts, numZShifts,
                                          tol, originalPoints, box, radii, baseStruc,
                                          randOrient)
            return list(filledAtom), "molecule"

    if isinstance(numMol, int) and numMol != 0:
        numXShifts = math.ceil(numMol ** (1/3))
        numYShifts = math.ceil(numMol / (numXShifts ** 2))
        numZShifts = math.ceil(numMol / (numXShifts * numYShifts))

        # Check if structure is monatomic or molecule
        if len(originalPoints) == 1:
            filledAtom = atomsDefiniteBox(numMol, maxattempts,
                                          originalPoints, box, radii, tol)
            return list(filledAtom), "atom"

        if len(originalPoints) > 1:
            filledAtom = moleculesDefiniteBox(numMol, maxattempts,
                                              originalPoints, box, radii, tol, baseStruc,
                                              randOrient)
            return list(filledAtom), "molecule"

    return "ERROR: No atoms found"

def drawMolSphere(struc, tol, radius, center, maxattempts, numMol, baseStruc, randOrient):
    '''Utilized to place molecules in a sphere'''
    sphere = Sphere3d(float(center[0]), float(center[1]), float(center[2]), radius)
    radii = setAtomicRadius()

    originalPoints = makeBase(struc)

    if not isinstance(numMol, int) and str(numMol).lower() != "fill":
        numMol = int(numMol)

    if str(numMol).lower() == 'fill':
        numShifts = math.ceil(2*sphere.radius / tol)

        # Check if structure is monatomic or molecule
        if len(originalPoints) == 1:
            filledAtom = atomFillSphere(numShifts, sphere,
                                        originalPoints, radii, tol)
            return list(filledAtom), "atom"

        if len(originalPoints) > 1:
            filledAtom = moleculeFillSphere(numShifts, sphere,
                                            originalPoints, radii, tol, baseStruc,
                                            randOrient)
            return list(filledAtom), "molecule"

    if isinstance(numMol, int) and numMol != 0:
        # Check if structure is monatomic or molecule
        if len(originalPoints) == 1:
            filledAtom = atomDefiniteSphere(numMol, maxattempts,
                                            originalPoints, sphere, radii, tol)
            return list(filledAtom), "atom"

        if len(originalPoints) > 1:
            filledAtom = moleculeDefiniteSphere(numMol, maxattempts,
                                                originalPoints, sphere, radii, tol, baseStruc,
                                                randOrient)
            return list(filledAtom), "molecule"

    return "ERROR: No atoms found"
