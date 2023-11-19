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
Currently only basic shapes such as boxes and spheres are generatable, but more complex meshes are in the works!

## Approach

TurtleMol takes an input molecules (or molecules) and copies them into a grid filling out a user defined volume.
The number of molecules in the volume can be specially defined, calculated via a desired density, or for filling space.

The orientation of the molecules can be randomized so structures are closer to their equilibrium point, but will still need outside optimization.
Molecules can also be placed randomly instead of in a grid, allowing for a more disordered structure.


## Current Features

- Generate a box or sphere of molecules
- Generate a box or sphere of molecules around an existing molecular structure
- Fill a volume, or place a specfied number of molecules in space
- Randomly orient molecules to better represent an equilibrium structure
- Read/Write structures from XYZ and PDB formats

## Features in Progress

- Generate structures based on meshes defined from outside software/packagaes such as Blender

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

## Citation

If you find this code helpful, please consider referencing it! We don't currently have a released article to cite,
but any reference to our work helps acknowledge the effort put into developing and maintaining this code base, 
provides support for further development!

## License

Distributed under the MIT License. See [LICENSE](https://github.com/Dfilono/TurtleMol/blob/main/LICENSE) for more information.

## Documentation

Documentation is in progress.

<img src="https://github.com/Dfilono/TurtleMol/blob/main/docs/images/logo.png">
