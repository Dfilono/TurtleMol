'''Tests for package'''

# pylin: skip-file

from TurtleMol.drawMol import drawMolBox, drawMolSphere
from TurtleMol.readWriteFiles import readStrucFile
from TurtleMol.defaultParams import defaultParams

def testReadWritePDB():
    iparamsSphere = {
        'shape' : 'sphere',
        'radius' : 5.0,
        'numMolecules' : 'fill',
        'structureFile' : "TurtleMol/examples/Water/water.pdb"
    }

    iparamsCube = {
        'shape' : 'cube',
        'sideLength' : 10.0,
        'numMolecules' : 'fill',
        'structureFile' : "TurtleMol/examples/Water/water.pdb"
    }

    dparams = defaultParams()

    for name in dparams:
        if name not in iparamsSphere:
            iparamsSphere[name] = dparams[name]

        if name not in iparamsCube:
            iparamsCube[name] = dparams[name]

    struc = readStrucFile(iparamsSphere['structureFile'])

    if iparamsSphere['baseStrucFile']:
        baseStruc = readStrucFile(iparamsSphere['baseStrucFile'])
    else:
        baseStruc = None

    outStrucCube, strucTypeCube = drawMolBox(struc, baseStruc, iparamsCube)
    assert len(outStrucCube) != 0, "Output structure should have some length"
    assert strucTypeCube == "molecule", "Structure type should be an atom"
    assert isinstance(outStrucCube[0][0][-1], str), 'Final place should be residue ID'

    outStrucSphere, strucTypeSphere = drawMolSphere(struc, baseStruc, iparamsSphere)
    assert len(outStrucSphere) != 0, "Output structure should have some length"
    assert strucTypeSphere == "molecule", "Structure type should be an atom"
    assert isinstance(outStrucSphere[0][0][-1], str), 'Final place should be residue ID'