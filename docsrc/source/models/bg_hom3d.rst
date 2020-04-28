.. highlight:: matlab
.. _bg_hom3d:

***********************
:mod:`bg_hom3d`
***********************

Background due to a homogeneous spin distribution in 3D


-----------------------------


Syntax
=========================================

.. code-block:: matlab

        info = bg_hom3d()
        P = bg_hom3d(r,param)
        P = bg_hom3d(r,param,lambda)

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


The expression for this model is :math:`B(t) = \mathrm{exp}\left(-\frac{8\pi^2}{9\sqrt{3}}\lambda c D |t|\right)` with :math:`D=\frac{\mu_0}{4\pi}\frac{(g_\mathrm{e}\mu_\mathrm{B})^2}{\hbar}`, with the concentration entered in spins/m\ :sup:`3`.

============= ============= ========= ============= ============= =============================================
 Variable       Symbol        Default   Lower bound   Upper bound      Description
============= ============= ========= ============= ============= =============================================
``param(1)``   :math:`c`       50          0.01          5000          Spin concentration (μM)
============= ============= ========= ============= ============= =============================================

-----------------------------


Description
=========================================

.. code-block:: matlab

        info = bg_hom3d()

Returns an ``info`` structure containing the specifics of the model:

* ``info.model`` -- Full name of the parametric model.
* ``info.nparam`` -- Total number of adjustable parameters.
* ``info.parameters`` -- Structure array with information on individual parameters.

-----------------------------


.. code-block:: matlab

    B = bg_hom3d(t,param)

Computes the background model ``B`` from the axis ``t`` according to the parameters array ``param`` for a modulation amplitude ``lambda=1``. The required parameters can also be found in the ``info`` structure.

-----------------------------

.. code-block:: matlab

    B = bg_hom3d(t,param,lambda)

Computes the background model ``B`` for a given modulation amplitude ``lambda`` (between 0 and 1).