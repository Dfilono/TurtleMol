import math
import sys
import pandas as pd
from Box3d import Box3d
from Sphere3D import Sphere3d
from shiftBox import AtomsFillBox, AtomsDefiniteBox, MoleculesFillBox, MoleculesDefiniteBox
from shiftSphere import AtomFillSphere, AtomDefiniteSphere, MoleculeDefiniteSphere, MoleculeFillSphere
from readWriteFiles import getElementData

'''
These functions are used to repeat a pattern of a given molecule or molecules. The goal is to repeat the pattern exactly the specified number of times
within the constraints of a the defined shape. A 'Fill' option was also added so that one could just define a space,
and fill said space with provided structure as many times as it can fit.
'''

def drawMolBox(struc, tol, dims, maxattempts, numMol):
    box = Box3d(0, 0, 0, dims)
    radii = setAtomicRadius()

    atoms = struc['Atom'].values.tolist()
    ogX = struc['X'].values.tolist()
    ogY = struc['Y'].values.tolist()
    ogZ = struc['Z'].values.tolist()

    original_points = [(atoms[i], ogX[i], ogY[i], ogZ[i]) for i in range(len(atoms))]

    if type(numMol) is not int and str(numMol).lower() != "fill":
        numMol = int(numMol)

    if str(numMol).lower() == "fill":
        num_x_shifts = math.ceil(box.length / tol)
        num_y_shifts = math.ceil(box.height / tol)
        num_z_shifts = math.ceil(box.width / tol)

        # Check if structure is monatomic or molecule
        if len(original_points) == 1:
            filledAtom = AtomsFillBox(num_x_shifts, num_y_shifts, num_z_shifts, tol, original_points, box, radii)
            return list(filledAtom), "atom"
        
        elif len(original_points) > 1:
            filledAtom = MoleculesFillBox(num_x_shifts, num_y_shifts, num_z_shifts, tol, original_points, box, radii)
            return list(filledAtom), "molecule"

        else:
            print('No atoms found')
            sys.exit()
        
        return list(filledAtom)

    elif type(numMol) is int and numMol != 0:
        num_x_shifts = math.ceil(numMol ** (1/3))
        num_y_shifts = math.ceil(numMol / (num_x_shifts ** 2))
        num_z_shifts = math.ceil(numMol / (num_x_shifts * num_y_shifts))

        # Check if structure is monatomic or molecule
        if len(original_points) == 1:
            filledAtom = AtomsDefiniteBox(numMol, maxattempts, original_points, box, radii, tol)
            return list(filledAtom), "atom"

        elif len(original_points) > 1:
            filledAtom = MoleculesDefiniteBox(numMol, maxattempts, original_points, box, radii, tol)
            return list(filledAtom), "molecule"

        else:
            print('No atoms found')
            sys.exit()

def drawMolSphere(struc, tol, radius, center, maxattempts, numMol):
    sphere = Sphere3d(float(center[0]), float(center[1]), float(center[2]), radius)
    radii = setAtomicRadius()

    atoms = struc['Atom'].values.tolist()
    ogX = struc['X'].values.tolist()
    ogY = struc['Y'].values.tolist()
    ogZ = struc['Z'].values.tolist()

    original_points = [(atoms[i], ogX[i], ogY[i], ogZ[i]) for i in range(len(atoms))]

    if type(numMol) is not int and str(numMol).lower() != "fill":
        numMol = int(numMol)

    if str(numMol).lower() == 'fill':
        num_shifts = math.ceil(2*sphere.radius / tol)

        # Check if structure is monatomic or molecule
        if len(original_points) == 1:
            filledAtom = AtomFillSphere(num_shifts, sphere, original_points, radii, tol)
            return list(filledAtom), "atom"

        elif len(original_points) > 1:
            filledAtom = MoleculeFillSphere(num_shifts, sphere, original_points, radii, tol)
            return list(filledAtom), "molecule"

        else:
            print('No atoms found')
            sys.exit()
    
    elif type(numMol) is int and numMol != 0:
        # Check if structure is monatomic or molecule
        if len(original_points) == 1:
            filledAtom = AtomDefiniteSphere(numMol, maxattempts, original_points, sphere, radii, tol)
            return list(filledAtom), "atom"

        elif len(original_points) > 1:
            filledAtom = MoleculeDefiniteSphere(numMol, maxattempts, original_points, sphere, radii, tol)
            return list(filledAtom), "molecule"

        else:
            print('No atoms found')
            sys.exit()

def setAtomicRadius():
    radii = {}

    ele = getElementData('AtomicRadius')
    
    radii = pd.Series(ele.AtomicRadius.values, index=ele.Symbol).to_dict()

    # Convert radii to Angstroms
    for k, v in radii.items():
        if v != 'nan':
            radii[k] = v * 0.01

    return radii





