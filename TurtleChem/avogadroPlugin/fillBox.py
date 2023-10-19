import argparse
import json
import sys
import tempfile
import pandas as pd
import TurtleChem

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

    if len(iparams) == 0:
        fd1, name1 = tempfile.mkstemp('_NoInput.txt')
        with open(name1, 'w') as f:
            for k,v in iparams:
                f.write(f'No input parameters')
            f.close()

    try:
        dparams = TurtleChem.defaultParams()
    except:
        fd1, name1 = tempfile.mkstemp("_ParamEror.txt")
        with open(name1, 'w') as f:
            f.write("defaultParams not found")
            f.close()

    for name in dparams:
        if name not in iparams:
            iparams[name] = dparams[name]

    return iparams

def runCommand():
    stdinStr = sys.stdin.read()
    opts = json.loads(stdinStr)
    iparams = generateParams(opts)

    if len(iparams) == 0:
        fd1, name1 = tempfile.mkstemp('_NoInput.txt')
        with open(name1, 'w') as f:
            for k,v in iparams:
                f.write(f'No input parameters')
            f.close()

    struc = pd.pd.read_csv(iparams['structureFile'], delim_whitespace=True,
                           skiprows=2, names=["Atom", "X", "Y", "Z"])
    
    try:
        if len(iparams['baseStrucFile']) > 0:
            baseStruc = TurtleChem.readStrucFile(iparams['baseStrucFile'])
        else:
            baseStruc = None
    except:
        fd1, name1 = tempfile.mkstemp("_ReadError.txt")
        with open(name1, 'w') as f:
            f.write("readWriteFiles not found")
            f.close()

    try:
        outStruc, strucType = TurtleChem.drawMolBox(struc, baseStruc, iparams)
    except:
        fd1, name1 = tempfile.mkstemp("_DrawError.txt")
        with open(name1, 'w') as f:
            f.write("drawMol not found")
            f.close()

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
