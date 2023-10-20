Input Parameters
================

The inputs for Turtle can be entered via the command line or an input file. Not all 
parameters need to entered, Turtle will rely on default values we have predefined. Some of those 
default values SHOULD be changed to ensure they meet the user's system, 
and those will be marked in the following list of parameters. The units for all lengths
are in Angstroms.

- shape: The shape of the object you would like to generate
    - Current options: cube, box, sphere
    - Default: box
- sideLength: The length of all sides in a cube
    - Default: 1.0
- XLen: The length of a non-cubical box (distance in the x-direction)
    - Default: 1.0
- YLen: The width of a non-cubical box (distance in the y-direction)
    - Default: 1.0
- ZLen: The height of a non-cubical box (distance in the z-direction)
    - Default: 1.0
- numMolecules: The number of molecules to place in the given volume
    - Current options: 
        - fill: A special parameter to fill the volume with maximum amount of molecules
        - Any integer
    - Default: fill
- tol: The minimum distance between two molecules in the space
    - Default: 1.0
- structureFile: The path to the XYZ or PDB file where your system is defined. This is REQUIRED
- baseStrucFile: The path to the XYZ or PDB file of an additional structure you want at the center of the new system
    - Default: None
- randomizeOrient: A boolean term that randomly orients the new molecules in the volume
    - Options: True, False
    - Default: False
- randFill: A boolean term that places molecules randomly in the volume rather than in a grid-like fashion. Only works when density or number of molecules is defined.
    - Options: True, False
    - Default: False
- density: Define the denisty of the system, and Turtle will calculate the number of molecules to place
    - Default: None
- atomRadius: Choose what type of radius you want the atoms to have
    - Options: AtomicRadius, CovalentRadius, VanDerWaalsRadius

Important parameters to keep in mind include denisty and atomRadius. Results may be unphysical for low denisty molecular systems 
if not specifically defined.

To build an input file, simply write a .txt file with the parameter=value. For example:

..code-block::
    shape=sphere
    radius=5.0
    tol=0.5
    randomizeOrient=True
    numMolecules=fill
    structureFile=path/to/file.txt
..code-block

Spaces between parameters are not allowed, spaces may be placed around the = sign. 

To run Turtle, open command line and input the following:

..code-block::
    python TurtleChem -i ./path/to/input/file.txt
..code-block
