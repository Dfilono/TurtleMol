Spherical volumes
=================

Use ``shape=sphere`` and ``radius`` to define the volume. The default center is
the origin; the CLI accepts another center through ``--center X Y Z``.

.. code-block:: text

   shape=sphere
   radius=15
   structureFile=water.xyz
   numMolecules=80
   randFill=True
   randomizeOrient=True
   outputFile=sphere.pdb

Atoms are tested against the spherical boundary with their selected atomic
radii. A requested count may be impossible when the radius is small or the
separation is large.
