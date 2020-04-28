.. highlight:: matlab
.. _bg_hom3d:

***********************
:mod:`bg_hom3d`
***********************

Multi-pulse DEER background in a homogenous medium


-----------------------------


Syntax
=========================================

.. code-block:: matlab

        info = bg_hom3d()
        P = bg_hom3d(r,param)
        P = bg_hom3d(r,param,lambda)

Parameters
    *   ``t`` - Time axis (N-array)
    *   ``param`` - Model parameters
    *   ``lambda`` -Dipolar pathway amplitude

Returns
    *   ``B`` - Model background (N-array)
    *   ``info`` - Model information (struct)


-----------------------------

Model
=========================================


============= ============= ========= ============= ============= ==================================
 Variable       Symbol        Default   Lower bound   Upper bound      Description
============= ============= ========= ============= ============= ==================================
``param(1)``   :math:`c`       50          0.01          5000          Spin concentration [`mol/dm^3`]
============= ============= ========= ============= ============= ==================================

-----------------------------


Description
=========================================

.. code-block:: matlab

        info = bg_hom3d()

Returns an ``info`` structure containing the specifics of the model:

* ``info.model`` -  Full name of the parametric model.
* ``info.nparam`` -  Total number of adjustable parameters.
* ``info.parameters`` - Structure array with information on individual parameters.

-----------------------------


.. code-block:: matlab

    B = bg_hom3d(t,param)

Computes the background model ``B`` from the axis ``t`` according to the parameters array ``param`` for a dipolar pathway amplitude ``lambda=1``. The required parameters can also be found in the ``info`` structure.

-----------------------------

.. code-block:: matlab

    B = bg_hom3d(t,param,lambda)

Computes the background model ``B`` for a given dipolar pathway amplitude ``lambda``.