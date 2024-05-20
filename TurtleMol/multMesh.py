'''Module handles multiple meshes/structures'''

import itertools
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
            
            if iparams['unitCells']:
                tol = 0
                iparams['unitCell'] = [iparams['unitCells'][i][0], iparams['unitCells'][i][1], iparams['unitCells'][i][2]]
                coord, strucType, cellParam = drawMolMesh(struc, baseStruc, iparams)
                cellParams.append(cellParam)
            else:
                tol = float(iparams['tol'])
                coord, strucType = drawMolMesh(struc, baseStruc, iparams)
                cellParams = None
            
            for mol in coord:
                newMol = []
                for atom in mol:
                    if len(atom) == 4:
                        atomData = (atom[0].capitalize(), atom[1], atom[2], atom[3])
                    if len(atom) == 5:
                        atomData = (atom[0].capitalize(), atom[1], atom[2], atom[3], atom[4])
                    newMol.append(atomData)
                
                if (kdTree is None or not isOverlapMoleculeKDTree(newMol, kdTree, indexToAtom, radii, tol)):
                    coords.append(newMol)

                    # Rebuild KDTree with newly added atoms
                    kdTree, indexToAtom = buildKDTreeMapping(coords, radii)
            
            strucTypes.append(strucType)

    elif isinstance(iparams['structureFile'], str):
            for i in range(len(meshList)):
                iparams['mesh'] = str(meshList[i])
                iparams['meshScale'] = float(scaleList[i])
                if iparams['unitCell']:
                    tol = 0
                    coord, strucType, cellParam = drawMolMesh(strucs, baseStruc, iparams)
                    cellParams = cellParam
                else:
                    tol = float(iparams['tol'])
                    coord, strucType = drawMolMesh(strucs, baseStruc, iparams)
                    cellParams = None

                for mol in coord:
                    newMol = []
                    for atom in mol:
                        if len(atom) == 4:
                            atomData = (atom[0].capitalize(), atom[1], atom[2], atom[3])
                        if len(atom) == 5:
                            atomData = (atom[0].capitalize(), atom[1], atom[2], atom[3], atom[4])
                        newMol.append(atomData)

                    if (kdTree is None or not isOverlapMoleculeKDTree(newMol, kdTree, indexToAtom, radii, tol)):
                        coords.append(newMol)

                        # Rebuild KDTree with newly added atoms
                        kdTree, indexToAtom = buildKDTreeMapping(coords, radii)
                strucTypes.append(strucType)

    return coords, 'molecule', cellParams