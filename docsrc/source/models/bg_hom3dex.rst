.. highlight:: matlab
.. _bg_hom3dex:

***********************
:mod:`bg_hom3dex`
***********************

Background due to a homogeneous spin distribution in 3D, with excluded volume

-----------------------------

Syntax
=========================================

.. code-block:: matlab

        info = bg_hom3dex()
        P = bg_hom3dex(r,param)
        P = bg_hom3dex(r,param,lambda)

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

This implements a hard-shell excluded-volume model, with distance of closest approach ``R`` (first parameter) and spin concentration ``c`` (second parameter, in μM).

The analytical expression for the decay function as a function of ``R`` and ``c`` is complicated, see `Kattnig et al, J.Phys.Chem.B 2013, 117, 16542 <https://pubs.acs.org/doi/abs/10.1021/jp408338q>`_.

============= =================== ========= ============= ============= ================================================
 Variable      Symbol              Default   Lower bound   Upper bound      Description
============= =================== ========= ============= ============= ================================================
``param(1)``    :math:`c`              50         0.01          1000          Spin concentration (μM)
``param(2)``    :math:`R`              1          0.1            20           Exclusion distance (nm)
============= =================== ========= ============= ============= ================================================

-----------------------------


Description
=========================================

.. code-block:: matlab

        info = bg_hom3dex()

Returns an ``info`` structure containing the specifics of the model:

* ``info.model`` -- Full name of the parametric model.
* ``info.nparam`` -- Total number of adjustable parameters.
* ``info.parameters`` -- Structure array with information on individual parameters.

-----------------------------


.. code-block:: matlab

    B = bg_hom3dex(t,param)

Computes the background model ``B`` from the axis ``t`` according to the parameters array ``param`` for a modulation amplitude ``lambda=1``. The required parameters can also be found in the ``info`` structure.

-----------------------------

.. code-block:: matlab

    B = bg_hom3dex(t,param,lambda)

Computes the background model ``B`` for a given modulation amplitude ``lambda`` (between 0 and 1).
