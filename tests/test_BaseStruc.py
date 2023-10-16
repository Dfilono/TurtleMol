'''Tests for package'''

# pylin: skip-file

from moleculeBuilder.drawMol import drawMolBox, drawMolSphere
from moleculeBuilder.readWriteFiles import readStrucFile
from moleculeBuilder.defaultParams import defaultParams

def testBaseStruc():
    iparamsSphere = {
        'shape' : 'sphere',
        'radius' : 5.0,
        'numMolecules' : 'fill',
        'structureFile' : "moleculeBuilder/examples/Methane_Water/water.xyz",
        'baseStrucFile' : "moleculeBuilder/examples/Methane_Water/ch4.xyz"
    }

    iparamsCube = {
        'shape' : 'cube',
        'sideLength' : 10.0,
        'numMolecules' : 'fill',
        'structureFile' : "moleculeBuilder/examples/Methane_Water/water.xyz",
        'baseStrucFile' : "moleculeBuilder/examples/Methane_Water/ch4.xyz"
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

    outStrucSphere, strucTypeSphere = drawMolSphere(struc, baseStruc, iparamsSphere)
    assert len(outStrucSphere) != 0, "Output structure should have some length"
    assert strucTypeSphere == "molecule", "Structure type should be an atom"
