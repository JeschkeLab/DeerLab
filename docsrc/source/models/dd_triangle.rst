.. highlight:: matlab
.. _dd_triangle:


***********************
:mod:`dd_triangle`
***********************

Triangle distribution

-----------------------------


Syntax
=========================================

.. code-block:: matlab

        info = dd_triangle()
        P = dd_triangle(r,param)

Parameters
    *   ``r`` - Distance axis (N-array)
    *   ``param`` - Model parameters
Returns
    *   ``P`` - Distance distribution (N-array)
    *   ``info`` - Model information (struct)

-----------------------------

Model
=========================================


This provides a simple triangular distribution.

============== ======================== ========= ============= ============= ========================
 Variable       Symbol                    Default   Lower bound   Upper bound      Description
============== ======================== ========= ============= ============= ========================
``param(1)``   :math:`r_0`                3.5       1.0              20         Maximum
``param(2)``   :math:`w_\mathrm{L}`       0.3       0.1              5          Left width
``param(3)``   :math:`w_\mathrm{R}`       0.3       0.1              5          Right width
============== ======================== ========= ============= ============= ========================


Example using default parameters:

.. image:: ../images/model_dd_triangle.png
   :width: 650px


-----------------------------


Description
=========================================

.. code-block:: matlab

        info = dd_triangle()

Returns an ``info`` structure containing the specifics of the model:

* ``info.model`` -  Full name of the parametric model.
* ``info.nparam`` -  Total number of adjustable parameters.
* ``info.parameters`` - Structure array with information on individual parameters.

-----------------------------


.. code-block:: matlab

    P = dd_triangle(r,param)

Computes the distance distribution model ``P`` from the axis ``r`` according to the parameters array ``param``. The required parameters can also be found in the ``info`` structure.

