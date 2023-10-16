'''
This module checks if atoms within a molecule overlaps
or if two different molecules overlap
'''

import math

def isOverlapMolecule(newAtom, filledAtoms, radii, tol):
    '''Checks if two molecules overlap'''
    for molecule in filledAtoms:
        for atom1, atom2 in zip(newAtom, molecule):
            distance = math.sqrt(
                (atom1[1] - atom2[1])**2 +
                (atom1[2] - atom2[2])**2 +
                (atom1[3] - atom2[3])**2
            )

            if distance < (radii[atom1[0]] + radii[atom2[0]]) + tol:
                return True
    return False

def isOverlapAtom(newAtom, filledAtoms, radii, tol):
    '''Checks if two atoms overlap'''
    for atom in filledAtoms:
        distance = math.sqrt(
            (newAtom[1] - atom[1])**2 +
            (newAtom[2] - atom[2])**2 +
            (newAtom[3] - atom[3])**2
        )

        if distance < (radii[newAtom[0]] + radii[atom[0]]) + tol:
            return True
    return False
