Troubleshooting
===============

No output file
--------------

Confirm that ``outputFile`` is set and ends in ``.xyz``, ``.pdb``, or
``.cjson``. Also verify that ``structureFile`` exists and has a supported
extension.

Too few molecules
-----------------

The requested count may not fit. Increase the volume or ``maxAttempts``, or
reduce ``tol``, the copy count, or the selected atomic radius. Random placement
becomes progressively less likely to succeed as the volume fills.

Unexpected command-line booleans
--------------------------------

Use an input file for boolean settings. With the current CLI parser,
``--randFill False`` is still a non-empty string and therefore evaluates true.

Mesh placement fails
--------------------

Check that the mesh is watertight, has sensible units and scale, and contains
no self-intersections. The ``rtree`` dependency must also be available for
trimesh spatial queries.

Unphysical contacts
-------------------

Try a larger ``tol`` or a larger radius model such as
``VanDerWaalsRadius``. TurtleMol does not perform force-field optimization, so
the result should be minimized before simulation.
