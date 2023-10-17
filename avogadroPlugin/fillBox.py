import argparse
import json
import sys
import pandas as pd
from moleculeBuilder import drawMol, \
    defaultParams, readWriteFiles

debug = True

def getOptions():
    boxOptions = {
        'Tolerance' : {
            'label' : 'Tolerance',
            'type' : 'float',
            'default' : 1.0,
            'precision' : 2,
            'toolTip' : 'Minimum distance between two molecules',
            'suffix' : 'Angstroms'
        },

        'Length' : {
            'label' : 'Length (X)',
            'type' : 'float',
            'default' : 1.0,
            'precision' : 2,
            'toolTip' : 'Length of Box',
            'suffix' : 'Angstroms'
        },

        'Width' : {
            'label' : 'Width (Y)',
            'type' : 'float',
            'default' : 1.0,
            'precision' : 2,
            'toolTip' : 'Width of Box',
            'suffix' : 'Angstroms'
        },

        'Height' : {
            'label' : 'Height (Z)',
            'type' : 'float',
            'default' : 1.0,
            'precision' : 2,
            'toolTip' : 'Height of Box',
            'suffix' : 'Angstroms'
        },

        'Randomize' : {
            'label' : 'Randomzie Orientation',
            'type' : 'boolean',
            'default' : False,
            'toolTip' : 'Randomize the orientation of molecules as they are placed'
        },

        'ogStruc' : {
            'label' : 'Original Structure',
            'type' : 'string',
            'default' : '',
            'toolTip' : 'Path to structure file for placement in center of box'
        }
    }

    opts = {'userOptions' : boxOptions}
    opts['inputMoleculeFormat'] = 'xyz'

    return opts

def generateParams(opts):
    iparams = {
        'shape' : 'box',
        'tol' : opts['Tolerance'],
        'Xlen' : opts['Length'],
        'Ylen' : opts['Width'],
        'Zlen' : opts['Height'],
        'numMolecules' : 'fill',
        'randomizeOrient' : opts['Randomize'],
        'structureFile' : opts['xyz'],
        'baseStrucFile' : opts['ogStruc']
    }

    dparams = defaultParams()
    for name in dparams:
        if name not in iparams:
            iparams[name] = dparams[name]

    return iparams

def runCommand():
    stdinStr = sys.stdin.read()
    opts = json.loads(stdinStr)
    iparams = generateParams(opts)

    struc = pd.pd.read_csv(iparams['structureFile'], delim_whitespace=True,
                           skiprows=2, names=["Atom", "X", "Y", "Z"])
    if len(iparams['baseStrucFile']) > 0:
        baseStruc = readWriteFiles.readStrucFile(iparams['baseStrucFile'])
    else:
        baseStruc = None

    outStruc, strucType = drawMol.drawMolBox(struc, baseStruc, iparams)

    columns = ['Atom', 'X', 'Y', 'Z']
    df = pd.DataFrame(columns=columns)

    if strucType == "molecule":
        for mol in outStruc:
            for atom in mol:
                df = df.append(pd.Series(atom, index=columns), ignore_index=True)
    else:
        df = pd.DataFrame(outStruc, columns=['Atom', 'X', 'Y', 'X'])

    file = f'{len(outStruc)}\n'

    for row in df:
        file += f'{row[0]}    {row[1]}    {row[2]}    {row[3]}\n'

    result = {}
    result['append'] = True
    result['moleculeFormat'] = 'xyz'
    result['xyz'] = file

    return result

if __name__ == "__main__":
    parser = argparse.ArgumentParser('TURTLE')
    parser.add_argument('--debug', action='store_true')
    parser.add_argument('--print-options', action='store_true')
    parser.add_argument('--run-command', action='store_true')
    parser.add_argument('--display-name', action='store_true')
    parser.add_argument('--menu-path', action='store_true')
    parser.add_argument('--lang', nargs='?', default='en')
    args = vars(parser.parse_args())

    debug = args['debug']

    if args['display_name']:
        print("Box")
    if args['menu_path']:
        print("&Extensions|TURTLE")
    if args['print_options']:
        print(json.dumps(getOptions()))
    elif args['run_command']:
        print(json.dumps(runCommand()))
