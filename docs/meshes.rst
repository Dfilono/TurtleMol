Mesh volumes
============

Set ``shape=mesh`` and provide a mesh path. TurtleMol loads the mesh through
trimesh, scales it, translates its minimum bound to the origin, and places
copies inside it.

.. code-block:: text

   shape=mesh
   mesh=container.stl
   meshScale=2.0
   structureFile=water.xyz
   numMolecules=200
   randFill=True
   maxAttempts=50000
   outputFile=mesh-fill.xyz

Use a closed, watertight mesh. Holes, self-intersections, inconsistent normals,
or unexpected source units can make containment tests unreliable. Inspect the
scaled mesh and output together when diagnosing missing or misplaced copies.

The Python parameter dictionary also supports ``onSurface`` and
``alignNormal`` for molecular tiles. Multi-mesh generation additionally needs
one tile, mesh, and global transformation matrix per component; it is an
advanced API workflow via :func:`TurtleMol.buildMultiMesh`.
