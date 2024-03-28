'''
This module checks if atoms within a molecule overlaps
or if two different molecules overlap
'''

import math
import numpy as np

def isOverlapMolecule(newAtom, filledAtoms, radii, tol): # Depreciated!!!!
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
    print('isOverlapMolecule is depreciated and will be removed from future versions!\n')
    return False

def isOverlapAtom(newAtom, filledAtoms, radii, tol): # Depreciated!!!!
    '''Checks if two atoms overlap'''
    for atom in filledAtoms:
        distance = math.sqrt(
            (newAtom[1] - atom[1])**2 +
            (newAtom[2] - atom[2])**2 +
            (newAtom[3] - atom[3])**2
        )

        if distance < (radii[newAtom[0]] + radii[atom[0]]) + tol:
            return True
        
    print('isOverlapMolecule is depreciated and will be removed from future versions!\n')
    return False

def isOverlapAtomKDTree(newAtom, kdTree, radii, tol):
    '''Checks if two atoms overlap using a KDTree'''
    atomName, x, y, z = newAtom[:4]
    radius = radii[str.capitalize(atomName)]
    nearbyAtoms = kdTree.query_ball_point([x, y, z], radius + tol)

    for i in nearbyAtoms:
        otherAtom = kdTree.data[i]
        distance = np.linalg.norm(np.array([x, y, z]) - otherAtom)

        if distance < (radius + radii[str.capitalize(kdTree.data[i][0])] + tol):
            return True
    return False

def isOverlapMoleculeKDTree(newMol, kdTree, radii, tol):
    '''Checks if any two atoms from a new molecule overlap using a KDTree'''
    for atom in newMol:
        atomName, x, y, z = atom[:4]
        radius = radii[str.capitalize(atomName)]
        nearbyAtoms = kdTree.query_ball_point([x, y, z], radius + tol)

        for i in nearbyAtoms:
            otherAtom = kdTree.data[i]
            distance = np.linalg.norm(np.array([x, y, z]) - otherAtom)

            if distance < (radius + radii[str.capitalize(kdTree.data[i][0])] + tol):
                return True
    return False