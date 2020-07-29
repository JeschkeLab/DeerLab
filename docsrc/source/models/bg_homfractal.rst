.. highlight:: matlab
.. _bg_homfractal:

***********************
:mod:`bg_homfractal`
***********************

Background due to a homogeneous spin distribution in a fractal dimension


-----------------------------


Syntax
=========================================

.. code-block:: matlab

        info = bg_homfractal()
        P = bg_homfractal(t,param)
        P = bg_homfractal(t,param,lambda)

Inputs
    *   ``t`` -- Time axis (N-array)
    *   ``param`` -- Model parameters
    *   ``lambda`` -- Modulation amplitude (between 0 and 1)

Outputs
    *   ``B`` -- Model background (N-array)
    *   ``info`` -- Model information (struct)


-----------------------------

Model
=========================================

This implements the background due to a homogeneous distribution of spins in a d-dimensional space, with d-dimensional spin concentration ``c_d``.


============= ============= ========= ============= ============= ===========================================
 Variable       Symbol        Default   Lower bound   Upper bound      Description
============= ============= ========= ============= ============= ===========================================
``param(1)``   :math:`c_d`     50          0.01          5000          Spin concentration (μmol/dm\ :sup:`d`)
``param(2)``   :math:`d`       3           0                6          Fractal dimension
============= ============= ========= ============= ============= ===========================================

-----------------------------


Description
=========================================

.. code-block:: matlab

        info = bg_homfractal()

Returns an ``info`` table containing the information of the model parameters and boundaries.

* ``info.model`` -- Full name of the parametric model.
* ``info.nparam`` -- Total number of adjustable parameters.
* ``info.parameters`` -- Structure array with information on individual parameters.

-----------------------------


.. code-block:: matlab

    B = bg_homfractal(t,param)

Computes the background model ``B`` from the axis ``t`` according to the parameters array ``param`` for a modulation amplitude ``lambda=1``. The required parameters can also be found in the ``info`` structure.

-----------------------------

.. code-block:: matlab

    B = bg_homfractal(t,param,lambda)

Computes the background model ``B`` for a given modulation amplitude ``lambda``.