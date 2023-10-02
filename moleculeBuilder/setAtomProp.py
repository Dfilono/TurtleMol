'''Sets properties of atoms'''

import pandas as pd
from readWriteFiles import getElementData

def setAtomicRadius():
    '''Defines the dataframe to find the atomic radius'''
    radii = {}

    ele = getElementData('AtomicRadius')
    
    radii = pd.Series(ele.AtomicRadius.values, index=ele.Symbol).to_dict()

    # Convert radii to Angstroms
    for k, v in radii.items():
        if v != 'nan':
            radii[k] = v * 0.01

    return radii

def setAtomicMass():
    '''Sets the mass of an atom'''
    mass = {}
    ele = getElementData('AtomicMass')
    mass = pd.Series(ele.AtomicMass.values, index=ele.Symbol).to_dict()

    return mass