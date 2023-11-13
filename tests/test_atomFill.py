'''Tests for package'''

# pylin: skip-file

from TurtleMol.drawMol import drawMolBox, drawMolSphere
from TurtleMol.readWriteFiles import readStrucFile
from TurtleMol.defaultParams import defaultParams

def testAtomFill():
    iparamsSphere = {
        'shape' : 'sphere',
        'radius' : 5.0,
        'numMolecules' : 'fill',
        'structureFile' : "TurtleMol/examples/Argon/argon.xyz"
    }

    iparamsCube = {
        'shape' : 'cube',
        'sideLength' : 10.0,
        'numMolecules' : 'fill',
        'structureFile' : "TurtleMol/examples/Argon/argon.xyz"
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

    print(iparamsCube)

    outStrucCube, strucTypeCube = drawMolBox(struc, baseStruc, iparamsCube)
    assert len(outStrucCube) != 0, "Output structure should have some length"
    assert strucTypeCube == "atom", "Structure type should be an atom"

    outStrucSphere, strucTypeSphere = drawMolSphere(struc, baseStruc, iparamsSphere)
    assert len(outStrucSphere) != 0, "Output structure should have some length"
    assert strucTypeSphere == "atom", "Structure type should be an atom"