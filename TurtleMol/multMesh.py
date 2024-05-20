'''Module handles multiple meshes/structures'''

import numpy as np
from .drawMol import drawMolMesh
from .isOverlap import isOverlapMoleculeKDTree, buildKDTreeMapping
from .setAtomProp import setAtomicRadius

def buildMultiMesh(strucs, baseStruc, iparams):
    '''Builds meshes of molecules based on a given number of meshes and atomic structures'''
    assert isinstance(iparams['mesh'], list), 'Only one mesh provided, or meshes not in list'
    assert isinstance(iparams['meshScale'], list), 'Please provide a scale for every provided mesh'
    
    coords = []
    strucTypes = []
    cellParams = []
    meshList = iparams['mesh']
    scaleList = iparams['meshScale']
    radii = setAtomicRadius(iparams['atomRadius'])

    # Create KD-tree for filledAtoms
    kdTree, indexToAtom = buildKDTreeMapping(coords, radii)

    if isinstance(iparams['structureFile'], list):
        assert len(strucs) == len(meshList), "If more than one structure, the number of structures needs to be the same as the number of meshes"

        for i in range(len(meshList)):
            iparams['mesh'] = str(meshList[i])
            iparams['meshScale'] = float(scaleList[i])
            struc = strucs[i]
            matrix = iparams['globalMatrix'][i]
            
            if iparams['unitCells']:
                tol = 0
                iparams['unitCell'] = [iparams['unitCells'][i][0], iparams['unitCells'][i][1], iparams['unitCells'][i][2]]
                coord, strucType, cellParam = drawMolMesh(struc, baseStruc, iparams)
                cellParams.append(cellParam)
            else:
                tol = float(iparams['tol'])
                coord, strucType = drawMolMesh(struc, baseStruc, iparams)
                cellParams = None

            allMolMesh = []
            for mol in coord:
                newMol = []
                for atom in mol:
                    if len(atom) == 4:
                        atomData = (atom[0].capitalize(), atom[1], atom[2], atom[3])
                    if len(atom) == 5:
                        atomData = (atom[0].capitalize(), atom[1], atom[2], atom[3], atom[4])
                    newMol.append(atomData)
                allMolMesh.extend(newMol)

            allMolMesh = applyGlobalTransform(allMolMesh, matrix)
            if (kdTree is None or not isOverlapMoleculeKDTree(allMolMesh, kdTree, indexToAtom, radii, tol)):
                coords.append(allMolMesh)
                # Rebuild KDTree with newly added atoms
                kdTree, indexToAtom = buildKDTreeMapping(coords, radii)
            
            strucTypes.append(strucType)

    elif isinstance(iparams['structureFile'], str):
            for i in range(len(meshList)):
                iparams['mesh'] = str(meshList[i])
                iparams['meshScale'] = float(scaleList[i])
                matrix = iparams['globalMatrix'][i]

                if iparams['unitCell']:
                    tol = 0
                    coord, strucType, cellParam = drawMolMesh(strucs, baseStruc, iparams)
                    cellParams = cellParam
                else:
                    tol = float(iparams['tol'])
                    coord, strucType = drawMolMesh(strucs, baseStruc, iparams)
                    cellParams = None

                allMolMesh = []
                for mol in coord:
                    newMol = []
                    for atom in mol:
                        if len(atom) == 4:
                            atomData = (atom[0].capitalize(), atom[1], atom[2], atom[3])
                        if len(atom) == 5:
                            atomData = (atom[0].capitalize(), atom[1], atom[2], atom[3], atom[4])
                        newMol.append(atomData)
                    allMolMesh.extend(newMol)

                allMolMesh = applyGlobalTransform(allMolMesh, matrix)
                if (kdTree is None or not isOverlapMoleculeKDTree(allMolMesh, kdTree, indexToAtom, radii, tol)):
                    coords.append(allMolMesh)

                    # Rebuild KDTree with newly added atoms
                    kdTree, indexToAtom = buildKDTreeMapping(coords, radii)
                strucTypes.append(strucType)

    return coords, 'molecule', cellParams


def applyGlobalTransform(data, matrix):
    '''Applies the global transformation to atom coordinates'''

    # Compute original centroid
    ogCentroid = computeCentroid(data)

    # Get Transformation Vector from Matrix
    transVector = matrix[:3, 3]

    # Compute Translation Vector
    transVector = transVector - np.array(ogCentroid)

    # Apply the translation to all atoms
    transformed = applyTranslation(data, transVector)
    
    return transformed

def computeCentroid(atoms):
    '''Computes the centroid of a set of atoms'''
    coords = np.array([atom[1:4] for atom in atoms])
    centroid = np.mean(coords, axis = 0)
    return centroid.tolist()

def applyTranslation(atoms, transVector):
    '''Applies the translation vector to all atoms'''
    translated = []
    for atom in atoms:
        translatedCoords = (np.array(atom[1:4]) + transVector).tolist()
        if len(atom) == 4:
            translated.append([atom[0]] + translatedCoords)
        elif len(atom) == 5:
            translated.append([atom[0]] + translatedCoords + [atom[4]])
    return translated
