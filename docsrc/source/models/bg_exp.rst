.. highlight:: matlab
.. _bg_exp:

***********************
:mod:`bg_exp`
***********************

Single-exponential background model

-----------------------------


Syntax
=========================================

.. code-block:: matlab

        info = bg_exp()
        P = bg_exp(t,param)
        P = bg_exp(t,param,lambda)

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

.. math::

   B(t) = \exp\left(-\lambda\kappa \vert t \vert\right)

============== =============== ========= ============= ============= ==============================
 Variable         Symbol        Default   Lower bound   Upper bound      Description
============== =============== ========= ============= ============= ==============================
``param(1)``   :math:`\kappa`   0.35         0            200          Decay rate
============== =============== ========= ============= ============= ==============================


Although the ``bg_exp`` model has the same functional form as ``bg_hom3d``, it is distinct since its parameter is a decay rate constant and not a spin concentration like for ``bg_hom3d``. 

-----------------------------


Description
=========================================

.. code-block:: matlab

        info = bg_exp()

Returns an ``info`` structure containing the specifics of the model:

* ``info.model`` -- Full name of the parametric model.
* ``info.nparam`` -- Total number of adjustable parameters.
* ``info.parameters`` -- Structure array with information on individual parameters.

-----------------------------


.. code-block:: matlab

    B = bg_exp(t,param)

Computes the background model ``B`` from the axis ``t`` according to the parameters array ``param`` for a modulation amplitude ``lambda=1``. The required parameters can also be found in the ``info`` structure.

-----------------------------

.. code-block:: matlab

    B = bg_exp(t,param,lambda)

Computes the background model ``B`` for a given modulation amplitude ``lambda``.

