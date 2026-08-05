'''Module handles multiple meshes/structures'''

import numpy as np
from .drawMol import drawMolMesh
from .setAtomProp import setAtomicRadius
from .makeStruc import applyGlobalTransform, applyTranslation, findMinPoint
from .treeManager import PendingKDManager

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
    rotAngleList = iparams['rotAngles']
    scaleList = []

    if isinstance(iparams['meshScale'], float):
        for i in meshList:
            scaleList.append(1)
    else:
        scaleList = iparams['meshScale']

    radii = setAtomicRadius(iparams['atomRadius'])

    # Create KD-tree for filledAtoms
    manager = PendingKDManager(radii, rebuildRate=500)

    def _matrixForIndex(i):
        '''Build per-mesh global matrix including optional per-axis scaling'''
        if scaleXList is not None and scaleYList is not None and scaleZList is not None:
            iparams['scaleX'] = float(scaleXList[i]) if scaleXList[i] is not None else scaleList[i]
            iparams['scaleY'] = float(scaleYList[i]) if scaleYList[i] is not None else scaleList[i]
            iparams['scaleZ'] = float(scaleZList[i]) if scaleZList[i] is not None else scaleList[i]

            scaleFactors = np.array([iparams['scaleX'], iparams['scaleY'], iparams['scalyZ']], dtype=float)

            scalingMatrix = np.eye(4, dtype=float)
            scalingMatrix[0,0] = scaleFactors[0]
            scalingMatrix[1,1] = scaleFactors[1]
            scalingMatrix[2,2] = scaleFactors[2]

            return np.dot(iparams['globalMatrix'][i], scalingMatrix)
        else:
            return iparams['globalMatrix'][i] * scaleList[i]
 
    def _flattenAndCapitalize(coord):
        '''Coord is list of molecules; flatten to list of atoms, capitalize types'''
        allAtoms = []
        for mol in coord:
            for atom in mol:
                if len(atom) == 4:
                    allAtoms.append((atom[0].capitalize(), atom[1], atom[2], atom[3]))
                elif len(atom) == 5:
                    allAtoms.append((atom[0].capitalize(), atom[1], atom[2], atom[3], atom[4]))
                else:
                    # Unexpected atom format; keep as-is but ensure type capitalization if possible
                    allAtoms.append(atom)
        return allAtoms
    
    def _tryAddStruc(allMolMesh, tol):
        """Overlap gate + commit to manager"""
        if tol is None or tol <= 0:
            coords.append(allMolMesh)
            manager.addMolecule(allMolMesh)
            manager.maybeRebuild()
            return True
        
        if not manager.overlapsMolecule(allMolMesh, tol):
            coords.append(allMolMesh)
            manager.addMolecule(allMolMesh)
            manager.maybeRebuild()
            return True
        
        return False

    if isinstance(iparams['structureFile'], list):
        assert len(strucs) == len(meshList), "If more than one structure, the number of structures needs to be the same as the number of meshes"

        for i in range(len(meshList)):
            iparams['mesh'] = str(meshList[i])
            iparams['meshScale'] = float(scaleList[i])
            iparams['rotAngles'] = rotAngleList[i]

            struc = strucs[i]
            matrix = _matrixForIndex(i)

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

            allMolMesh = _flattenAndCapitalize(coord)
            allMolMesh = applyGlobalTransform(all, matrix)

            _tryAddStruc(allMolMesh, tol)
            strucTypes.append(strucType)

    elif isinstance(iparams['structureFile'], str):
            for i in range(len(meshList)):
                iparams['mesh'] = str(meshList[i])
                iparams['meshScale'] = float(scaleList[i])
                iparams['rotAngles'] = rotAngleList[i]

                struc = strucs
                matrix = _matrixForIndex(i)

                if iparams['unitCell']:
                    tol = 0
                    coord, strucType, cellParam = drawMolMesh(strucs, baseStruc, iparams)
                    cellParams.append(cellParam)
                else:
                    tol = float(iparams['tol'])
                    coord, strucType = drawMolMesh(strucs, baseStruc, iparams)
                    cellParams = None

                allMolMesh = _flattenAndCapitalize(coord)
                allMolMesh = applyGlobalTransform(allMolMesh, matrix)

                _tryAddStruc(allMolMesh, tol)
                strucTypes.append(strucType)

    transVector = np.array([0, 0, 0]) - findMinPoint(coords)
    coords = applyTranslation(coords, transVector, 'molecule')

    return coords, 'molecule', cellParams

