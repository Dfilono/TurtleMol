Boxes and cubes
===============

Use ``shape=box`` with ``Xlen``, ``Ylen``, and ``Zlen`` for a rectangular
volume. Use ``shape=cube`` with ``sideLength`` for equal edges.

.. code-block:: text

   shape=box
   Xlen=30
   Ylen=20
   Zlen=15
   structureFile=water.xyz
   numMolecules=fill
   tol=1.0
   outputFile=box.xyz

The box begins at the origin. Grid filling places as many non-overlapping
copies as allowed by the volume and tolerance. Set a finite ``numMolecules``
and ``randFill=True`` for random positions; increase ``maxAttempts`` if a dense
job stops before reaching its requested count.
