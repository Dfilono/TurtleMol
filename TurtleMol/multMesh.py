'''Module handles multiple meshes/structures'''

import numpy as np
from .drawMol import drawMolMesh
from .isOverlap import isOverlapMoleculeKDTree, buildKDTreeMapping
from .setAtomProp import setAtomicRadius

def buildMultiMesh(strucs, baseStruc, iparams):
    '''Builds meshes of molecules based on a given number of meshes and atomic structures'''
    assert isinstance(iparams['mesh'], list), 'Only one mesh provided, or meshes not in list'
    #assert isinstance(iparams['meshScale'], list), 'Please provide a scale for every provided mesh'
    
    coords = []
    strucTypes = []
    cellParams = []
    meshList = iparams['mesh']
    scaleXList = iparams['scaleX']
    scaleYList = iparams['scaleY']
    scaleZList = iparams['scaleZ']
    scaleList = []

    if isinstance(iparams['meshScale'], float):
        for i in meshList:
            scaleList.append(1)
    else:
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
            if scaleXList is not None and scaleYList is not None and scaleZList is not None:
                iparams['scaleX'] = float(scaleXList[i]) if scaleXList[i] is not None else scaleList[i]
                iparams['scaleY'] = float(scaleYList[i]) if scaleYList[i] is not None else scaleList[i]
                iparams['scaleZ'] = float(scaleZList[i]) if scaleZList[i] is not None else scaleList[i]
            
                scaleFactors = np.array([iparams['scaleX'], iparams['scaleY'], iparams['scaleZ']])

                # Create scaling matrix
                scalingMatrix = np.eye(4)
                scalingMatrix[0, 0] = scaleFactors[0]
                scalingMatrix[1, 1] = scaleFactors[1]
                scalingMatrix[2, 2] = scaleFactors[2]

                matrix = np.dot(iparams['globalMatrix'][i], scalingMatrix)
            else:
                matrix = iparams['globalMatrix'][i] * scaleList[i]
            
            if iparams['unitCells'][i] is not None:
                tol = 0
                iparams['unitCell'] = [iparams['unitCells'][i][0], iparams['unitCells'][i][1], iparams['unitCells'][i][2]]
                iparams['angle'] = [iparams['angles'][i][0], iparams['angles'][i][1], iparams['angles'][i][2]]
                coord, strucType, cellParam = drawMolMesh(struc, baseStruc, iparams)
                cellParams.append(cellParam)
            else:
                tol = float(iparams['tol'])
                iparams['unitCell'] = None
                coord, strucType = drawMolMesh(struc, baseStruc, iparams)
                cellParams.append(None)

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
                struc = strucs[i]
                if scaleXList is not None and scaleYList is not None and scaleZList is not None:
                    iparams['scaleX'] = float(scaleXList[i]) if scaleXList[i] is not None else scaleList[i]
                    iparams['scaleY'] = float(scaleYList[i]) if scaleYList[i] is not None else scaleList[i]
                    iparams['scaleZ'] = float(scaleZList[i]) if scaleZList[i] is not None else scaleList[i]

                    scaleFactors = np.array([iparams['scaleX'], iparams['scaleY'], iparams['scaleZ']])

                    # Create scaling matrix
                    scalingMatrix = np.eye(4)
                    scalingMatrix[0, 0] = scaleFactors[0]
                    scalingMatrix[1, 1] = scaleFactors[1]
                    scalingMatrix[2, 2] = scaleFactors[2]

                    matrix = np.dot(iparams['globalMatrix'][i], scalingMatrix)
                else:
                    matrix = iparams['globalMatrix'][i] * scaleList[i]

                if iparams['unitCell']:
                    tol = 0
                    coord, strucType, cellParam = drawMolMesh(strucs, baseStruc, iparams)
                    cellParams.append(cellParam)
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

    transVector = np.array([0, 0, 0]) - findMinPoint(coords)
    coords = applyTranslation(coords, transVector, 'molecule')
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

def applyTranslation(atoms, transVector, strucType='atom'):
    '''Applies the translation vector to all atoms'''
    translated = []
    if strucType == 'atom':
        for atom in atoms:
            translatedCoords = (np.array(atom[1:4]) + transVector).tolist()
            if len(atom) == 4:
                translated.append([atom[0]] + translatedCoords)
            elif len(atom) == 5:
                translated.append([atom[0]] + translatedCoords + [atom[4]])
    elif strucType == 'molecule':
        for mol in atoms:
            newMol = []
            for atom in mol:
                translatedCoords = (np.array(atom[1:4]) + transVector).tolist()
                if len(atom) == 4:
                    newMol.append([atom[0]] + translatedCoords)
                elif len(atom) == 5:
                    newMol.append([atom[0]] + translatedCoords + [atom[4]])
            translated.append(newMol)
    return translated

def findMinPoint(mols):
    '''Finds the minimum point of the total structure'''
    coords = np.array([atom[1:4] for mol in mols for atom in mol])
    minPoint = np.min(coords, axis=0)
    return minPoint