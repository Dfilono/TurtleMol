'''Tests for package'''

# pylin: skip-file

from moleculeBuilder.drawMol import atomRandSphere, atomsRandBox
from moleculeBuilder.readWriteFiles import readStrucFile
from moleculeBuilder.defaultParams import defaultParams

def testAtomFill():
    iparamsSphere = {
        'shape' : 'sphere',
        'radius' : 5.0,
        'numMolecules' : '20',
        'randFill' : True,
        'structureFile' : "../examples/Argon/argon.xyz"
    }

    iparamsCube = {
        'shape' : 'cube',
        'sideLength' : 10.0,
        'numMolecules' : '20',
        'randFill' : True,
        'structureFile' : "../examples/Argon/argon.xyz"
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

    outStrucCube, strucTypeCube = atomsRandBox(struc, baseStruc, iparamsCube)
    assert len(outStrucCube) != 0, "Output structure should have some length"
    assert strucTypeCube == "atom", "Structure type should be an atom"

    outStrucSphere, strucTypeSphere = atomRandSphere(struc, baseStruc, iparamsSphere)
    assert len(outStrucSphere) != 0, "Output structure should have some length"
    assert strucTypeSphere == "atom", "Structure type should be an atom"