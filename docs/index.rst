TurtleMol documentation
=======================

TurtleMol builds molecular starting structures by repeating an atom, molecule,
or crystallographic unit cell inside a box, sphere, or mesh. It can fill the
available volume, place a requested number of copies, or estimate the number of
copies from a target mass density.

TurtleMol creates *starting geometries*. Always inspect and energy-minimize
generated structures before using them in a simulation.

Quick example
-------------

.. code-block:: console

   python -m pip install TurtleMol
   TurtleMol --structureFile water.xyz --shape cube --sideLength 20 \
             --numMolecules 100 --tol 1.0 --outputFile water-box.xyz

The ``TurtleMol`` command and ``python -m TurtleMol`` are equivalent.

.. toctree::
   :maxdepth: 2
   :caption: User guide

   getting-started
   command-line
   input
   file-formats
   box
   sphere
   meshes
   rho
   out
   troubleshooting

.. toctree::
   :maxdepth: 2
   :caption: Python API

   python-api

Project links
-------------

* `Source code <https://github.com/Dfilono/TurtleMol>`_
* `Issue tracker <https://github.com/Dfilono/TurtleMol/issues>`_
* `MIT license <https://github.com/Dfilono/TurtleMol/blob/main/LICENSE>`_
