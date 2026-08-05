Output
======

Set ``outputFile`` to a path ending in ``.xyz``, ``.pdb``, or ``.cjson``.
TurtleMol writes no structure file when the option is omitted.

XYZ is the simplest and most portable output but contains only elements and
coordinates. PDB can include crystallographic cell and connectivity records in
unit-cell workflows. CJSON stores atomic numbers and a flattened coordinate
array, with optional unit-cell values.

Generated structures may contain fewer copies than requested when boundary and
overlap constraints make placement impossible. Always validate the output atom
count and geometry before passing it to downstream software.
