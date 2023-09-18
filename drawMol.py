import math
import pandas as pd
from Box3d import Box3d
from readWriteFiles import getElementData

'''
This function is used to repeat a pattern of a given molecule or molecules. The goal is to repeat the pattern exactly the specified number of times
within the constraints of a the defined shape. A 'Fill' option was also added so that one could just define a space,
and fill said space with provided structure as many times as it can fit.
'''

def drawMol(struc, tol, dims, fromWall, numMol): # TODO add a way to account for atomic radius of atoms
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

        filledAtom = shiftAtomsFill(num_x_shifts, num_y_shifts, num_z_shifts, tol, original_points, box, radii)
        
        return list(filledAtom)

    elif type(numMol) is int and numMol != 0:

        filledAtom = shiftAtomsDefinite(numMol, tol, original_points, box, radii)
        
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

def is_overlap(new_atom, filled_atoms, radii):
    for atom in filled_atoms:
        distance = math.sqrt(
            (new_atom[1] - atom[1])**2 +
            (new_atom[2] - atom[2])**2 +
            (new_atom[3] - atom[3])**2
        )

        if distance < (radii[new_atom[0]] + radii[atom[0]]):
            return True
    return False

def shiftAtomsFill(x, y, z, tol, og, box, radii):
    filledAtom = []
    filledPositions = set()
    for zshift in range(z):
        for yshift in range(y):
            for xshift in range(x):
                for atom in og:
                    # Calculate the shift for each tile and new point
                    new_x = atom[1] + xshift * tol
                    new_y = atom[2] + yshift * tol
                    new_z = atom[3] + zshift * tol

                    # Adjust for atomic radiss
                    new_x_min = new_x - radii[atom[0]]
                    new_y_min = new_y - radii[atom[0]]
                    new_z_min = new_z - radii[atom[0]]
                    new_x_max = new_x + radii[atom[0]]
                    new_y_max = new_y + radii[atom[0]]
                    new_z_max = new_z + radii[atom[0]]

                    # Check if the new atom fits within the box
                    if  (
                        new_x_min >= box.x and new_x_max <= box.x + box.length and
                        new_y_min >= box.y and new_y_max <= box.y + box.width and
                        new_z_min >= box.z and new_z_max <= box.z + box.height
                    ):
                        new_atom = (atom[0], new_x, new_y, new_z)

                if not is_overlap(new_atom, filledAtom, radii):
                    filledAtom.append(new_atom)
                    filledPositions.add(new_atom)

    return filledAtom

def shiftAtomsDefinite(numMol, tol, og, box, radii): # NOTE Currently does not work as intended
    filledAtom = []
    filledPositions = set()

    for i in range(numMol):
        for j in range(numMol):
            for k in range(numMol):
                for atom in og:
                    # Calculate the shift for each tile and new point
                    new_x = atom[1] + k*tol
                    new_y = atom[2] + j*tol
                    new_z = atom[3] + i*tol

                    # Adjust for atomic radiss
                    new_x_min = new_x - radii[atom[0]]
                    new_y_min = new_y - radii[atom[0]]
                    new_z_min = new_z - radii[atom[0]]
                    new_x_max = new_x + radii[atom[0]]
                    new_y_max = new_y + radii[atom[0]]
                    new_z_max = new_z + radii[atom[0]]

                    # Check if the new atom fits within the box
                    if  (
                        new_x_min >= box.x and new_x_max <= box.x + box.length and
                        new_y_min >= box.y and new_y_max <= box.y + box.width and
                        new_z_min >= box.z and new_z_max <= box.z + box.height
                    ):
                        new_atom = (atom[0], new_x, new_y, new_z)

                if not is_overlap(new_atom, filledAtom, radii):
                    filledAtom.append(new_atom)
                    filledPositions.add(new_atom)

    return filledAtom