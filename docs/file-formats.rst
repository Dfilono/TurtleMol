File formats
============

TurtleMol selects a format from the path extension (case-insensitive).

XYZ
---

XYZ input must have an atom count, a comment line, then whitespace-separated
``element x y z`` records. Coordinates are interpreted as ångströms. XYZ does
not preserve bonds, residues, or unit-cell metadata.

PDB
---

TurtleMol reads ``ATOM``, ``HETATM``, ``CRYST1``, and ``CONECT`` records. It
writes generated atoms as ``HETATM`` records and can preserve cell information
for unit-cell workflows. Element symbols should be present in columns 77-78;
the atom name is used as a fallback.

CJSON
-----

Chemical JSON input uses ``atoms.elements.number`` and flattened
``atoms.coords.3d`` arrays. An optional ``unitCell`` object may contain ``a``,
``b``, ``c``, ``alpha``, ``beta``, and ``gamma``.

Meshes
------

Mesh files are loaded by trimesh. Common formats include STL, OBJ, and PLY.
Closed, watertight meshes give the most reliable inside/outside tests. See
:doc:`meshes` for mesh-specific controls.
