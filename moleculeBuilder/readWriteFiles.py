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
    '''Reads structure file'''
    # Reads XYZ
    if filePath[-3:] == 'xyz':
        return pd.read_csv(filePath, delim_whitespace=True,
                           skiprows=2, names=["Atom", "X", "Y", "Z"])
    # Reads PDB
    if filePath[-3:] == 'pdb':
        return readPdb(filePath)
    
    return f"ERROR: Issue generating file to {filePath}"

def readPdb(filePath):
    '''Reads structure file if a pdb'''
    atoms = []
    residues = []

    with open(filePath, 'r') as pdbFile:
        for line in pdbFile:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                # Parse relevant fields from the PDB format
                atomName = line[12:16].strip()
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                residue = line[17:20].strip()

                atoms.append([atomName, x, y, z])
                residues.append(residue)

    df = pd.DataFrame(atoms, columns=['Atom', 'X', 'Y', 'Z'])
    df['Residue'] = residues

    return df

def writeOutput(data, filePath, strucType):
    '''Writes data to output file'''
    # Writes XYZ
    if filePath[-3:] == 'xyz':
        writeXYZ(data, filePath, strucType)
    # Writes PDB
    if filePath[-3:] == 'pdb':
        writePdb(data, filePath)
    else:
        return f"ERROR: Filetype {filePath[-3:]} not supported"

def writePdb(data, filePath):
    '''Writes a pdb file from results'''
    with open(filePath, 'w', encoding = 'utf-8') as pdbFile:
        molNum = 0
        resNum = 1
        for mol in data:
            for atom in mol:
                # Extract information for each atom
                atomName = atom[0]
                x = atom[1]
                y = atom[2]
                z = atom[3]
                residue = atom[4]

                # Format the line in PDB format
                line = f"ATOM  {(atomName + str(molNum)):<4} {(residue + str(molNum)):<3}   " \
                + str(resNum) + f"    {x:8.3f}{y:8.3f}{z:8.3f}\n"
                pdbFile.write(line)
            molNum += 1
            resNum += 1
        pdbFile.write("TER\n")

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
    Mass: g/mol
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
