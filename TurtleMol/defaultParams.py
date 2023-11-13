def defaultParams():
    '''
Keywords:
    shape = "Cube" or "Box" are the currently supported shapes
    sideLength = float(length of sides in a cube)
    Xlen = float(length of a box)
    Ylen = float(width of a box)
    Zlen = float(depth of a box)
    numMolecules = int(Number of molecules in box)
    tol = float(minimum distance between molecules)
    inputFile = Path to the input file if you are using one
    structureFile = Path to the structure file being used. 
        Required for usage. XYZ and PDB formats are currently supported
    outputFile = Path to output file if desired
    baseStrucFile = Path to file for a structure that is present, but
        not deuplicated
    baseStrucCenter = coordinates to place the geometric
        center of the base structure
    randomizeOrient = Boolean to randomize orientation of
        molecules
    randFill = When placing a set number of molecules, place them randomly
    desnity = total mass of molecules over the total volume
        of the system in g/mL
    atomRadius = Choose whether to use the AtomicRadius, CovalentRadius,
        or VanDerWaalsRadius

'''
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
        'maxAttempts' : 10000,
        'baseStrucCenter' : [0, 0, 0],
        'baseStrucFile' : None,
        'randomizeOrient' : False,
        'randFill' : False,
        'density' : None,
        'atomRadius' : 'AtomicRadius',
    }

    return params
