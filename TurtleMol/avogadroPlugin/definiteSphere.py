import argparse
import json
import sys
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import pandas as pd
import TurtleMol
import io

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

        'Radius' : {
            'label' : 'Radius',
            'type' : 'float',
            'default' : 1.0,
            'precision' : 2,
            'toolTip' : 'Radius of Sphere',
            'suffix' : 'Angstroms'
        },

        'Randomize' : {
            'label' : 'Randomzie Orientation',
            'type' : 'boolean',
            'default' : False,
            'toolTip' : 'Randomize the orientation of molecules as they are placed'
        },

        'ogStruc' : {
            'label' : 'Structure From File',
            'type' : 'filePath',
            'default' : '',
            'toolTip' : 'Path to structure file for placement in center of box'
        },

        'numMol' : {
            'label' : 'Number of Molecules',
            'type' : 'integer',
            'default' : 1,
            'toolTip' : 'Choose the number of molecules you want to add to the system!',
        },

        'Placement' : {
            'label' : 'Randomly Place',
            'type' : 'boolean',
            'default' : False,
            'toolTip' : 'Randomly place molecules instead of placing them in a grid'
        }
    }

    opts = {'userOptions' : boxOptions}
    opts['inputMoleculeFormat'] = 'xyz'

    return opts

def generateParams(opts):
    iparams = {
        'shape' : 'box',
        'tol' : opts['Tolerance'],
        'radius' : opts['Radius'],
        'numMolecules' : opts['numMol'],
        'randomizeOrient' : opts['Randomize'],
        'structureFile' : opts['xyz'],
        'baseStrucFile' : opts['ogStruc'],
        'randFill' : opts['Placement']
    }

    dparams = TurtleMol.defaultParams()

    for name in dparams:
        if name not in iparams:
            iparams[name] = dparams[name]

    return iparams

def runCommand():
    stdinStr = sys.stdin.read()
    opts = json.loads(stdinStr)
    iparams = generateParams(opts)

    strucIO = io.StringIO(iparams['structureFile'])

    struc = pd.read_csv(strucIO, delim_whitespace=',',
                            skiprows=2, names=["Atom", "X", "Y", "Z"])
    

    if len(iparams['baseStrucFile']) > 0:
        baseStruc = TurtleMol.readStrucFile(iparams['baseStrucFile'])
    else:
        baseStruc = None


    outStruc, strucType = TurtleMol.drawMolSphere(struc, baseStruc, iparams)

    columns = ['Atom', 'X', 'Y', 'Z']
    df = pd.DataFrame(columns=columns)

    if strucType == "molecule":
        for mol in outStruc:
            for atom in mol:
                df = pd.concat([df, pd.DataFrame([atom], columns=columns)], ignore_index=True)
    else:
        df = pd.DataFrame(outStruc, columns=['Atom', 'X', 'Y', 'X'])

    file = f'{len(df)}\n\n'

    if len(df) != 0:
        for index, row in df.iterrows():
            file += f"{row['Atom']:<2}  {row['X']:> 10.5f}   {row['Y']:> 8.5f}   {row['Z']:> 8.5f}\n"
    else:
        file = iparams['structureFile']

    result = {}
    result['append'] = False
    result['moleculeFormat'] = 'xyz'
    result['xyz'] = file

    return result

if __name__ == "__main__":
    parser = argparse.ArgumentParser('Choose Atoms in Sphere')
    parser.add_argument('--debug', action='store_true')
    parser.add_argument('--print-options', action='store_true')
    parser.add_argument('--run-command', action='store_true')
    parser.add_argument('--display-name', action='store_true')
    parser.add_argument('--menu-path', action='store_true')
    parser.add_argument('--lang', nargs='?', default='en')
    args = vars(parser.parse_args())

    debug = args['debug']

    if args['display_name']:
        print("Choose Number of Atoms Sphere")
    if args['menu_path']:
        print("&Extensions|TurtleMol")
    if args['print_options']:
        print(json.dumps(getOptions()))
    elif args['run_command']:
        print(json.dumps(runCommand()))
