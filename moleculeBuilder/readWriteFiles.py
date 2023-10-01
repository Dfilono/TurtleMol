'''
This module is for reading and writing files as needed. 
Since both the input and output files are optional, these functions can be avoided.
The structure file is required however, and currently the only supported format is xyz files.
Data on elements is also imported and searchable. Primary use at this time is for the atomic radius
'''

import os
import pandas as pd

def getInput(filePath):
    '''Reads input file if provided'''
    params = {}
    with open(filePath, encoding = 'utf-8') as f:
        for line in f:
            key, value = line.strip().split('=')
            params[key.strip()] = value.strip()

    return params

def readStrucFile(filePath):
    '''Reads structure file if an xyz'''
    return pd.read_csv(filePath, delim_whitespace=True,
                       skiprows=2, names=["Atom", "X", "Y", "Z"])

def writeXYZ(data, filePath, strucType):
    '''Writes an xyz file from results'''
    if strucType == 'molecule':
        columns = ['Atom', 'X', 'Y', 'Z']
        df = pd.DataFrame(columns=['Atom', 'X', 'Y', 'Z'])

        for mol in data:
            for atom in mol:
                df = df._append(pd.Series(atom, index=columns),
                                ignore_index = True)

        with open(filePath, 'w', encoding = 'utf-8') as f:
            f.write(str(len(df['Atom'])))
            f.write('\n\n')
            f.write(df.to_string(header=False, index=False))

    else:
        df = pd.DataFrame(data, columns=['Atom', 'X', 'Y', 'X'])
        with open(filePath, 'w', encoding = 'utf-8') as f:
            f.write(str(len(df['Atom'])))
            f.write('\n\n')
            f.write(df.to_string(header=False, index=False))

def getElementData(param):
    '''
    Import Atomic Data
    Database Information
    Denisty Units: g/mL, Temp Units: Kelvin, Heat Capacity Units: J/(g*K),
    Heat of Vaporization/Fusion Units = kJ/mol
    Molar Volume Units: mL, Thermal Conductivty Units: Watts/(m * K),
    Atomic/Covalent/Van Det Waals Radius Units = pm
    Ionization Energy Units = kJ/mol, Electron Affinity Units = kJ/mol
    Data From: https://github.com/Bluegrams/periodic-table-data
    '''

    dirname = os.path.dirname(__file__)
    filename = os.path.join(dirname, 'data/ElementData.csv')

    elements = pd.read_csv(filename, delimiter=",", header=0)
    columnsToRemove = ["DiscoveredBy", "Discovery", "AbundanceCrust", "AbundanceUniverse"]
    elements = elements.drop(labels=columnsToRemove, axis = 1)

    return pd.DataFrame(elements, columns=['Symbol', param])
