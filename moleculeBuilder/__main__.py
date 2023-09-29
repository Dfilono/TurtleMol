import sys
import argparse
from Box3d import drawBox
from drawMol import drawMolBox, drawMolSphere
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
        'sphereCenter' : [0, 0, 0],
        'radius' : 1.0,
        'numMolecules' : 1,
        'tol' : 1.0,
        'fromWall' : 1.0,
        'maxAttempts' : 10000
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
    parser.add_argument('-r', '--radius', type=float, help="Radius of a spehere in Angstroms", default=dparams['radius'])
    parser.add_argument('-center', '--center', nargs='+', type=float, help="X, Y, Z coordinates of the center of a sphere", default=dparams['sphereCenter'])
    parser.add_argument('-n', '--numMolecules', type=int, help="Number of molecules", default=dparams['numMolecules'])
    parser.add_argument('-t', '--tol', type=float, help="Minimum distance between molecules", default=dparams['tol'])
    parser.add_argument('-w', '--fromWall', type=float, help="Minimum distance from wall", default=dparams['fromWall'])
    parser.add_argument('-a', '--maxAttempts', type=int, help="Maximum iterations for finite sized systems", default=dparams['maxAttempts'])
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

    # Get structure
    struc = readStrucFile(iparams['structureFile'])
    print(struc)

    if iparams['shape'].lower() == 'box' or iparams['shape'].lower() == 'cube':
        # Define the box
        dims = drawBox(iparams)
        print(dims)

        # Generate the new structure
        out_struc = drawMolBox(struc, float(iparams['tol']), dims, float(iparams['maxAttempts']), iparams['numMolecules'])
        print(len(out_struc))
    
    elif iparams['shape'].lower() == 'sphere':
        # Define Sphere
        center = iparams['sphereCenter']
        radius = iparams['radius']

        # Generate the new structure
        out_struc = drawMolSphere(struc, float(iparams['tol']), float(radius), center, float(iparams['maxAttempts']), iparams['numMolecules'])

    if iparams['outputFile']:
        writeXYZ(out_struc, iparams['outputFile'])


if __name__ == "__main__":
    main()

