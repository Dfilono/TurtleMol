Getting started
===============

Requirements
------------

TurtleMol requires Python 3.7 or newer. Its runtime dependencies include
NumPy, pandas, SciPy, trimesh, and rtree.

Installation
------------

Install the released package from PyPI:

.. code-block:: console

   python -m pip install TurtleMol

To use the latest source checkout:

.. code-block:: console

   git clone https://github.com/Dfilono/TurtleMol.git
   cd TurtleMol
   python -m pip install .

Create a source molecule
------------------------

The input structure is the tile TurtleMol repeats. Save this as ``water.xyz``:

.. code-block:: text

   3
   water
   O  0.000000  0.000000  0.000000
   H  0.957200  0.000000  0.000000
   H -0.239987  0.927297  0.000000

Generate a structure
--------------------

.. code-block:: console

   TurtleMol --structureFile water.xyz --shape cube --sideLength 20 \
             --numMolecules 100 --tol 1.0 --outputFile water-box.xyz

The file extension selects XYZ, PDB, or CJSON format. Open the result in a
molecular viewer and check its boundaries and close contacts. The generated
structure is not geometry-optimized.
