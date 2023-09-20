import math
import sys
import pandas as pd
from Box3d import Box3d
from shiftAtoms import shiftAtomsFill, shiftAtomsDefinite
from shiftMols import shiftMoleculesFill, shiftMoleculesDefinite
from readWriteFiles import getElementData

'''
This function is used to repeat a pattern of a given molecule or molecules. The goal is to repeat the pattern exactly the specified number of times
within the constraints of a the defined shape. A 'Fill' option was also added so that one could just define a space,
and fill said space with provided structure as many times as it can fit.
'''

def drawMol(struc, tol, dims, maxattempts, numMol):
    box = Box3d(0, 0, 0, dims[0], dims[1], dims[2])
    radii = setAtomicRadius()

    atoms = struc['Atom'].values.tolist()
    ogX = struc['X'].values.tolist()
    ogY = struc['X'].values.tolist()
    ogZ = struc['X'].values.tolist()

    original_points = [(atoms[i], ogX[i] + radii[atoms[i]], ogY[i] + radii[atoms[i]], ogZ[i] + radii[atoms[i]]) for i in range(len(atoms))]

    if type(numMol) is not int and str(numMol).lower() != "fill":
        numMol = int(numMol)

    if str(numMol).lower() == "fill":
        num_x_shifts = math.ceil(box.length / tol)
        num_y_shifts = math.ceil(box.height / tol)
        num_z_shifts = math.ceil(box.width / tol)

        # Check if structure is monatomic or molecule
        if len(original_points[0]) == 1:
            filledAtom = shiftAtomsFill(num_x_shifts, num_y_shifts, num_z_shifts, tol, original_points, box, radii)

        elif len(original_points[0]) > 1:
            filledAtom = shiftMoleculesFill(num_x_shifts, num_y_shifts, num_z_shifts, tol, original_points, box, radii)

        else:
            print('No atoms found')
            sys.exit()
        
        return list(filledAtom)

    elif type(numMol) is int and numMol != 0:
        num_x_shifts = math.ceil(numMol ** (1/3))
        num_y_shifts = math.ceil(numMol / (num_x_shifts ** 2))
        num_z_shifts = math.ceil(numMol / (num_x_shifts * num_y_shifts))

        # Check if structure is monatomic or molecule
        if len(original_points[0]) == 1:
            filledAtom = shiftAtomsDefinite(numMol, maxattempts, original_points, box, radii)

        elif len(original_points[0]) > 1:
            filledAtom = shiftMoleculesDefinite(numMol, maxattempts, original_points, box, radii)

        else:
            print('No atoms found')
            sys.exit()
        
        return list(filledAtom)
                        
def setAtomicRadius():
    radii = {}

    ele = getElementData('AtomicRadius')
    
    radii = pd.Series(ele.AtomicRadius.values, index=ele.Symbol).to_dict()

    # Convert radii to Angstroms
    for k, v in radii.items():
        if v != 'nan':
            radii[k] = v * 0.01

    return radii





