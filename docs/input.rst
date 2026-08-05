Input files and parameters
==========================

An input file is UTF-8 text containing one ``key=value`` pair per line. Key
names are case-sensitive. Do not add comments or blank lines: every line is
parsed as an assignment.

.. code-block:: text

   shape=sphere
   radius=15.0
   tol=0.5
   randomizeOrient=True
   numMolecules=fill
   structureFile=water.xyz
   outputFile=water-sphere.pdb

Run it with ``TurtleMol --inputFile run.txt``.

Parameter reference
-------------------

.. list-table::
   :header-rows: 1
   :widths: 22 18 60

   * - Key
     - Default
     - Description
   * - ``shape``
     - ``box``
     - ``box``, ``cube``, ``sphere``, ``mesh``, or ``multimesh``.
   * - ``structureFile`` / ``outputFile``
     - required / none
     - Input tile and output structure paths.
   * - ``sideLength``
     - ``1.0``
     - Cube edge length in Å.
   * - ``Xlen``, ``Ylen``, ``Zlen``
     - ``1.0``
     - Box dimensions in Å.
   * - ``radius`` / ``sphereCenter``
     - ``1.0`` / ``[0,0,0]``
     - Sphere radius and center. Set list-valued options through Python; the
       text input parser does not deserialize comma-separated lists.
   * - ``numMolecules``
     - ``1``
     - Integer copy count or ``fill``.
   * - ``tol``
     - ``1.0``
     - Minimum separation for overlap checks, in Å.
   * - ``density``
     - none
     - Target density in g/mL. Overrides count-based placement.
   * - ``randomizeOrient`` / ``randFill``
     - ``False``
     - Randomize molecular orientations / positions.
   * - ``maxAttempts``
     - ``10000``
     - Attempt limit for random placement.
   * - ``atomRadius``
     - ``AtomicRadius``
     - Element radius used in overlap checks.
   * - ``baseStrucFile`` / ``baseStrucCenter``
     - none / ``[0,0,0]``
     - Non-repeated structure and its desired center.
   * - ``mesh`` / ``meshScale``
     - none / ``1.0``
     - Mesh path and uniform scale.
   * - ``unitCell`` / ``angle``
     - none / ``[90,90,90]``
     - Cell lengths and angles, normally read from input.
   * - ``hexagonal`` / ``rotAngles``
     - ``False`` / ``[0,0,0]``
     - Hexagonal packing and fixed tile rotation.
   * - ``onSurface`` / ``alignNormal``
     - ``False``
     - Mesh surface placement and face-normal alignment.

``scaleX``, ``scaleY``, ``scaleZ``, and ``padding`` are experimental.
``globalMatrixPath`` supplies transforms for multi-mesh generation.
