import sys
import argparse
from Box3d import drawBox
from drawMol import drawMol
from readWriteFiles import *

'''
Keywords:
    shape = "Cube" or "Box" are the currently supported shapes
    sideLength = float(length of sides in a cube)
    Xlen = float(length of a box)
    Ylen = float(width of a box)
    Zlen = float(depth of a box)
    numMolecules = int(Number of molecules in box)
    tol = float(minimum distance between molecules)
    fromWall = minimum distance from wall
    inputFile = Path to the input file if you are using one
    structureFile = Path to the structure file being used. Required for usage. Only the xyz format is currently supported
    outputFile = Path to output file if desired

'''

def defaultParams():
    params = {
        'shape' : 'box',
        'sideLength' : 1.0,
        'Xlen' : 1.0,
        'Ylen' : 1.0,
        'Zlen' : 1.0,
        'numMolecules' : 1,
        'tol' : 1.0,
        'fromWall' : 1.0,
    }

    return params

def parseCommandLine(dparams):
    parser = argparse.ArgumentParser(description='Tile Fill')

    # Command line arguemnnts
    parser.add_argument('-i', '--inputFile', type=str, help="Path to input file")
    parser.add_argument('-struc', '--structureFile', type=str, help='Path to structure file')
    parser.add_argument('-s', '--shape', type=str, help="Shape (Box or Cube)", default=dparams['shape'])
    parser.add_argument('-sl', '--sideLength', type=float, help="Dimensions of a cube in Angstroms", default=dparams['sideLength'])
    parser.add_argument('-xl', '--Xlen', type=float, help="X dimension of a box in Angstroms", default=dparams['Xlen'])
    parser.add_argument('-yl', '--Ylen', type=float, help="Y dimension of a box in Angstroms", default=dparams['Ylen'])
    parser.add_argument('-zl', '--Zlen', type=float, help="Z dimension of a box in Angstroms", default=dparams['Zlen'])
    parser.add_argument('-n', '--numMolecules', type=int, help="Number of molecules", default=dparams['numMolecules'])
    parser.add_argument('-t', '--tol', type=float, help="Minimum distance between molecules", default=dparams['tol'])
    parser.add_argument('-w', '--fromWall', type=float, help="Minimum distance from wall", default=dparams['fromWall'])
    parser.add_argument('-out', '--outputFile', type=str, help="Path for output file if desired")

    return parser.parse_args()

def main():
    # Init default params
    dparams = defaultParams()

    # Init arguement parser
    args = parseCommandLine(dparams)

    iparams = {}
    if args.inputFile:
        iparams = getInput(args.inputFile)
    else:
        for arg_name in vars(args):
            iparams[arg_name] = getattr(args, arg_name)

    for name in dparams:
        if name not in iparams:
            iparams[name] = dparams[name]

    if not iparams['structureFile']:
        print('Please provide an initial structure')
        sys.exit()
    
    print(iparams)

    # Define the box
    dims = drawBox(iparams)
    print(dims)

    # Get structure
    struc = readStrucFile(iparams['structureFile'])
    print(struc)

    # Generate the new structure
    out_struc = drawMol(struc, float(iparams['tol']), dims, float(iparams['fromWall']), iparams['numMolecules'])
    print(len(out_struc))

    if iparams['outputFile']:
        writeXYZ(out_struc, iparams['outputFile'])


if __name__ == "__main__":
    main()

