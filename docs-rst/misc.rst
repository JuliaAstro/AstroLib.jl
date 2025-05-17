cirrange
~~~~~~~~

.. function:: cirrange(number[, max]) -> restricted_number

Purpose
'''''''

Force a number into a given range :math:`[0, \text{max})`.

Argument
''''''''

-  ``number``: the number to modify. Can be a scalar or an array.
-  ``max`` (optional numerical argument): specify the extremum of the
   range :math:`[0, \text{max})` into which the number should be
   restricted. If omitted, defaults to ``360.0``.

Output
''''''

The converted number or array of numbers, as ``AbstractFloat``.

Example
'''''''

Restrict an array of numbers in the range :math:`[0, 2\pi)` as if they
are angles expressed in radians:

.. code-block:: julia

    cirrange([4pi, 10, -5.23], 2.0*pi)
    # => 3-element Array{Float64,1}:
    #     0.0
    #     3.71681
    #     1.05319

Notes
'''''

This function does not support the ``radians`` keyword like IDL
implementation. Use ``2.0*pi`` as second argument to restrict a number
to the same interval.

Code of this function is based on IDL Astronomy User's Library.

--------------

rad2sec
~~~~~~~

.. function:: rad2sec(rad) -> seconds

Purpose
'''''''

Convert from radians to seconds.

Argument
''''''''

-  ``rad``: number of radians. It can be either a scalar or a vector.

Output
''''''

The number of seconds corresponding to ``rad``. If ``rad`` is an array,
an array of the same length is returned.

Example
'''''''

.. code-block:: julia

    rad2sec(1)
    # => 206264.80624709636

Notes
'''''

Use ``sec2rad`` to convert seconds to radians.

--------------

sec2rad
~~~~~~~

.. function:: sec2rad(sec) -> radians

Purpose
'''''''

Convert from seconds to radians.

Argument
''''''''

-  ``sec``: number of seconds. It can be either a scalar or a vector.

Output
''''''

The number of radians corresponding to ``sec``. If ``sec`` is an array,
an array of the same length is returned.

Example
'''''''

.. code-block:: julia

    sec2rad(3600*30)
    # => 0.5235987755982988

Notes
'''''

Use ``rad2sec`` to convert radians to seconds.
