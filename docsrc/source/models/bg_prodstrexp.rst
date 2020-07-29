.. highlight:: matlab
.. _bg_prodstrexp:


***********************
:mod:`bg_prodstrexp`
***********************

Product of two stretched exponentials background model

-----------------------------


Syntax
=========================================

.. code-block:: matlab

        info = bg_prodstrexp()
        P = bg_prodstrexp(t,param)
        P = bg_prodstrexp(t,param,lambda)

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

:math:`B(t) = \exp\left(-\lambda\kappa_1 \vert t \vert^{d_1}\right) \exp\left(-\lambda\kappa_2 \vert t\vert^{d_2}\right)`

============= ================== ========= ============= ============= ==============================
 Variable         Symbol          Default   Lower bound   Upper bound      Description
============= ================== ========= ============= ============= ==============================
``param(1)``   :math:`\kappa_1`    3.5         0            200         1st strexp decay rate
``param(2)``   :math:`d_1`         1           0            6           1st strexp stretch factor
``param(3)``   :math:`\kappa_2`    3.5         0            200         2nd strexp decay rate
``param(4)``   :math:`d_2`         1           0            6           2nd strexp stretch factor
============= ================== ========= ============= ============= ==============================

-----------------------------


Description
=========================================

.. code-block:: matlab

        info = bg_prodstrexp()

Returns an ``info`` table containing the information of the model parameters and boundaries.

* ``info.model`` -- Full name of the parametric model.
* ``info.nparam`` -- Total number of adjustable parameters.
* ``info.parameters`` -- Structure array with information on individual parameters.

-----------------------------


.. code-block:: matlab

    B = bg_prodstrexp(t,param)

Computes the background model ``B`` from the axis ``t`` according to the parameters array ``param`` for a modulation amplitude ``lambda=1``. The required parameters can also be found in the ``info`` structure.

-----------------------------

.. code-block:: matlab

    B = bg_prodstrexp(t,param,lambda)

Computes the background model ``B`` for a given modulation amplitude ``lambda``.
