import pandas as pd

'''
This module is for reading and writing files as needed. Since both the input and output files are optional, these functions can be avoided.
The structure file is required however, and currently the only supported format is xyz files.
Data on elements is also imported and searchable. Primary use at this time is for the atomic radius
'''

def getInput(file_path):
    params = {}
    with open(file_path) as f:
        for line in f:
            key, value = line.strip().split('=')
            params[key.strip()] = value.strip()
    
    return params

def readStrucFile(file_path):
    return pd.read_csv(file_path, delim_whitespace=True, skiprows=2, names=["Atom", "X", "Y", "Z"])

def writeXYZ(data, filepath):
    columns = ['Atom', 'X', 'Y', 'Z']
    df = pd.DataFrame(columns=['Atom', 'X', 'Y', 'Z'])

    for mol in data:
        for atom in mol:
            df = df.append(pd.Series(atom, index=columns), ignore_index = True)

    with open(filepath, 'w') as f:
        f.write(str(len(df['Atom'])))
        f.write('\n\n')
        f.write(df.to_string(header=False, index=False))

def getElementData(param):
    # Import Atomic Data
    # Database Information
    # Denisty Units: g/mL, Temp Units: Kelvin, Heat Capacity Units: J/(g*K), Heat of Vaporization/Fusion Units = kJ/mol
    # Molar Volume Units: mL, Thermal Conductivty Units: Watts/(m * K), Atomic/Covalent/Van Det Waals Radius Units = pm
    # Ionization Energy Units = kJ/mol, Electron Affinity Units = kJ/mol
    # Data From: https://github.com/Bluegrams/periodic-table-data

    elements = pd.read_csv("data/ElementData.csv", delimiter=",", header=0)
    columnsToRemove = ["DiscoveredBy", "Discovery", "AbundanceCrust", "AbundanceUniverse"]
    elements = elements.drop(labels=columnsToRemove, axis = 1)

    return pd.DataFrame(elements, columns=['Symbol', param])
