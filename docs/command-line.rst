Command-line interface
======================

Run TurtleMol through its installed command or as a module:

.. code-block:: console

   TurtleMol [options]
   python -m TurtleMol [options]

Use ``TurtleMol --help`` for the options supported by the installed version.

Core options
------------

.. list-table::
   :header-rows: 1
   :widths: 27 18 55

   * - Option
     - Default
     - Meaning
   * - ``--structureFile``
     - required
     - XYZ, PDB, or CJSON tile to repeat.
   * - ``--outputFile``
     - none
     - Destination XYZ, PDB, or CJSON file.
   * - ``--shape``
     - ``box``
     - ``box``, ``cube``, ``sphere``, ``mesh``, or ``multimesh``.
   * - ``--numMolecules``
     - ``1``
     - Number of copies. Input files also accept ``fill``.
   * - ``--tol``
     - ``1.0``
     - Minimum separation in ångströms.
   * - ``--density``
     - none
     - Target density in g/mL; takes precedence over copy count.
   * - ``--randomizeOrient``
     - false
     - Randomize the orientation of molecular copies.
   * - ``--randFill``
     - false
     - Place copies randomly instead of on a grid.
   * - ``--maxAttempts``
     - ``10000``
     - Attempt limit for random filling.
   * - ``--atomRadius``
     - ``AtomicRadius``
     - ``AtomicRadius``, ``CovalentRadius``, or ``VanDerWaalsRadius``.
   * - ``--baseStrucFile``
     - none
     - Structure retained but not repeated.

Shape options
-------------

``box`` uses ``--Xlen``, ``--Ylen``, and ``--Zlen``. ``cube`` uses
``--sideLength``. ``sphere`` uses ``--radius`` and ``--center X Y Z``.
``mesh`` uses ``--mesh`` and optionally ``--meshScale``.

.. note::

   For reliable boolean values, especially false values, prefer an input file.
   The current CLI applies Python's ``bool`` conversion, under which any
   non-empty string is true.

When ``--inputFile`` is supplied, other CLI values are not merged into it.
Missing keys receive defaults. Density takes precedence over ``numMolecules``.
