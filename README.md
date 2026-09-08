<h1 align='center'>TurtleMol</h1>
<h4 align='center'>Design, create, and generate unique, complex molecular structures of any shape and size!</h4>



<p align="center">
    <a href="https://github.com/Dfilono/TurtleMol/actions/workflows/python-package.yml">
        <img src="https://github.com/Dfilono/TurtleMol/actions/workflows/python-package.yml/badge.svg" alt="Build Status ">
    </a>
    <a href="https://codecov.io/gh/Dfilono/TurtleMol">
        <img src="https://codecov.io/gh/Dfilono/TurtleMol/branch/main/graph/badge.svg?token=P643JEUWZC" alt="Codecov">
    </a>
    <a href="https://github.com/Dfilono/TurtleMol/blob/main/LICENSE" target="_blank">
        <img src="https://img.shields.io/github/license/Dfilono/TurtleMol" alt="License">
    </a>
    <a href="https://github.com/Dfilono/TurtleMol" target="_blank">
        <img src="https://img.shields.io/github/repo-size/Dfilono/TurtleMol" alt="Repo size">
    </a>
    <a href="https://github.com/psf/black" target="_blank">
        <img src="https://img.shields.io/badge/code%20style-black-000000.svg" alt="Black style">
    </a>
    <a href="https://github.com/PyCQA/pylint" target="_blank">
        <img src="https://img.shields.io/badge/linting-pylint-yellowgreen" alt="Black style">
    </a>
</p>

## Introduction

The goal of TurtleMol is to ease the process of
generating unique and complex molecular structures for computational chemistry. Available as a standalone python-package,
or as a plugin for [Avogadro2](https://www.openchemistry.org/projects/avogadro2/) TurtleMol aims to quickly generate
the intial state of molecular strcutres that are ready for optimization via quamtum or classical calculations!

Generating intial structure files for calculations can be a tedious process, one that TurtleMol aims to fix. 
Large structures can be generated quickly, and are very tunable based on the need of the user! 

## Approach

TurtleMol takes an input molecules (or molecules) and copies them into a grid filling out a user defined volume.
The number of molecules in the volume can be specially defined, calculated via a desired density, or for filling space.

The orientation of the molecules can be randomized so structures are closer to their equilibrium point, but will still need outside optimization.
Molecules can also be placed randomly instead of in a grid, allowing for a more disordered structure.


## Current Features

- Generate a box, sphere, or mesh of molecules
- Generate a box, sphere, or mesh of molecules around an existing molecular structure
- Fill a volume, or place a specfied number of molecules in space
- Randomly orient molecules to better represent an equilibrium structure
- Orient molecules to align with nearest normal surface vectors of a mesh
- Place molecules along the surface of a mesh
- Make supercell structure using unit cell parameters a, b, c, alpha, beta, gamma
- Read/Write structures from XYZ, PDB, cjson formats

## Installation

You can install the latest development version of TurtleMol from the [Github Repository](https://github.com/Dfilono/TurtleMol).

    git clone https://github.com/Dfilono/TurtleMol
    cd TurtleMol
    pip install .

Or you can download it from PyPi:

    pip install TurtleMol

TurtleMol is also available as a plugin for [Avogadro2](https://www.openchemistry.org/projects/avogadro2/)
and can be installed via the following instructions.

<img src="https://github.com/Dfilono/TurtleMol/blob/main/docs/images/installationPart1.png">

<img src="https://github.com/Dfilono/TurtleMol/blob/main/docs/images/installationPart2.png">

Note that for the plugin to function, the TurtleMol python package and its dependecies must also be installed
in the same Python environment that is referenced by Avogadro.

## Contributing

We have a lot of wishlist features that can be seen [HERE](https://github.com/Dfilono/TurtleMol/blob/main/WISHLIST.md). If you want to help add any of these features, or others you think TurtleMol would benefit from, let me know! Submit a pull request with your update (please test it in your own fork first), and I'll check it!

## Citation

If you find this code helpful, please consider citing it!

```
@article{
doi:10.26434/chemrxiv-2025-r4806,
author = {Dominick Filonowich  and Geoffrey Hutchison  and Christopher Wilmer },
title = {TurtleMol: Flexible Generation of Complex Molecular Systems for Computational Chemistry},
journal = {ChemRxiv},
volume = {2025},
number = {1128},
pages = {},
year = {2025},
doi = {10.26434/chemrxiv-2025-r4806},
URL = {https://chemrxiv.org/doi/abs/10.26434/chemrxiv-2025-r4806},
eprint = {https://chemrxiv.org/doi/pdf/10.26434/chemrxiv-2025-r4806},
abstract = {TurtleMol is an open-source Python package that aims to help users generate large, complex molec- ular systems. In the current version, users can generate systems by filling volumes defined by basic geometric shapes (e.g. cube, sphere), or by shapes of arbitrary gemoetries defined meshes created in other software (such as Blender, SOLIDWORKS, AutoDeskInventor). Volumes can be filled by user- defined patterns of atoms (e.g., a water molecule or unit-cell of quartz) that tile the specified volume. Several options are available for filling a system: (1) the tiles at fixed spacing, (2) the tiles at fixed density, and (3) tiles positioned/oriented randomly. TurtleMol does not optimize atomic positions in any way. The package is freely available on GitHub at https://github.com/Dfilono/TurtleMol}}
```

## License

Distributed under the MIT License. See [LICENSE](https://github.com/Dfilono/TurtleMol/blob/main/LICENSE) for more information.

## Documentation

The complete user guide and Python API reference live in [`docs/`](docs/index.rst).
Build the HTML site locally with:

    python -m sphinx -W -b html docs docs/_build/html

<img src="https://github.com/Dfilono/TurtleMol/blob/main/docs/images/logo.png">
