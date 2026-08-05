Density-based generation
========================

Set ``density`` in g/mL to let TurtleMol calculate a copy count from the tile's
molar mass and the target volume. Density mode takes precedence over
``numMolecules``.

.. code-block:: text

   shape=cube
   sideLength=30
   density=1.0
   structureFile=water.xyz
   randomizeOrient=True
   outputFile=water-density.xyz

The calculation depends on recognized element symbols and the selected
geometric volume. The requested bulk density does not guarantee a relaxed or
physically realistic local structure; minimize and equilibrate the output with
appropriate simulation software.
