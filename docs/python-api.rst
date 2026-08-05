Python API
==========

TurtleMol exposes its primary geometry, generation, file, and transformation
helpers at package level. Begin with :func:`TurtleMol.defaultParams`, override
the values needed by your job, and pass the complete dictionary to a draw
function.

.. code-block:: python

   import TurtleMol

   params = TurtleMol.defaultParams()
   params.update({
       'shape': 'cube',
       'sideLength': 20.0,
       'numMolecules': 100,
       'tol': 1.0,
       'randomizeOrient': True,
   })

   tile, unit_cell, connectivity = TurtleMol.readStrucFile('water.xyz')
   if unit_cell:
       params['unitCell'] = [unit_cell['a'], unit_cell['b'], unit_cell['c']]
       params['angle'] = [unit_cell['alpha'], unit_cell['beta'],
                          unit_cell['gamma']]

   result = TurtleMol.drawMolBox(tile, None, params,
                                 conect=connectivity)
   molecules, structure_type, *metadata = result
   TurtleMol.writeOutput(molecules, 'water-box.xyz', structure_type)

Generation
----------

.. py:function:: defaultParams()

   Return a new dictionary containing every supported generation parameter.
   The dictionary can be modified and passed to a draw function.
.. autofunction:: TurtleMol.drawMolBox
.. autofunction:: TurtleMol.drawMolSphere
.. autofunction:: TurtleMol.drawMolMesh
.. autofunction:: TurtleMol.buildMultiMesh

File handling
-------------

.. autofunction:: TurtleMol.getInput
.. autofunction:: TurtleMol.readStrucFile
.. autofunction:: TurtleMol.readPdb
.. autofunction:: TurtleMol.readMesh
.. autofunction:: TurtleMol.writeOutput
.. autofunction:: TurtleMol.writePdb
.. autofunction:: TurtleMol.writeXYZ
.. autofunction:: TurtleMol.getElementData

Geometry
--------

.. autoclass:: TurtleMol.Box3d
   :members:

.. autoclass:: TurtleMol.Sphere3d
   :members:

.. autoclass:: TurtleMol.mesh3D
   :members:

Transformations and calculations
--------------------------------

.. autofunction:: TurtleMol.calcCenter
.. autofunction:: TurtleMol.reCenter
.. autofunction:: TurtleMol.shiftPoints
.. autofunction:: TurtleMol.Reorient
.. autofunction:: TurtleMol.calcDensity
.. autofunction:: TurtleMol.calcNumMol
.. autofunction:: TurtleMol.applyGlobalTransform
.. autofunction:: TurtleMol.computeCentroid
.. autofunction:: TurtleMol.applyTranslation
.. autofunction:: TurtleMol.alignToNormal
.. autofunction:: TurtleMol.alignVectors
.. autofunction:: TurtleMol.placeOnSurfaceNormal

Return values
-------------

Draw functions normally return ``(structures, structure_type)``. Unit-cell
paths append cell metadata, and box generation with connectivity may append an
updated connectivity mapping. Code that supports both modes can unpack the
first two values and collect the remainder, as in the example above.
