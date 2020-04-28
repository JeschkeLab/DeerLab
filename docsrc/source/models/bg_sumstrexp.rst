.. highlight:: matlab
.. _bg_sumstrexp:


***********************
:mod:`bg_sumstrexp`
***********************

Sum of two stretched exponentials background parametric model

-----------------------------


Syntax
=========================================

.. code-block:: matlab

        info = bg_sumstrexp()
        P = bg_sumstrexp(r,param)
        P = bg_sumstrexp(r,param,lambda)

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

:math:`B(t) = A_1\exp \left(-\lambda\kappa_1 \vert t \vert^{d_1/3}\right) + (1-A_1)\exp\left(-\lambda\kappa_2 \vert t \vert^{d_2/3}\right)`

============= ================= ========= ============= ============= ==============================
 Variable       Symbol           Default   Lower bound   Upper bound      Description
============= ================= ========= ============= ============= ==============================
``param(1)``  :math:`\kappa_1`     3.5         0            200         1st strexp decay rate
``param(2)``  :math:`d_1`          3           0            6           1st strexp fractal dimension
``param(3)``  :math:`\kappa_2`     3.5         0            200         2nd strexp decay rate
``param(4)``  :math:`d_2`          3           0            6           2nd strexp fractal dimension
``param(5)``  :math:`A_1`          0.5         0            1           Relative amplitude
============= ================= ========= ============= ============= ==============================

-----------------------------


Description
=========================================

.. code-block:: matlab

        info = bg_sumstrexp()

Returns an ``info`` structure containing the specifics of the model:

* ``info.model`` -  Full name of the parametric model.
* ``info.nparam`` -  Total number of adjustable parameters.
* ``info.parameters`` - Structure array with information on individual parameters.

-----------------------------


.. code-block:: matlab

    B = bg_sumstrexp(t,param)

Computes the background model ``B`` from the axis ``t`` according to the parameters array ``param`` for a dipolar pathway amplitude ``lambda=1``. The required parameters can also be found in the ``info`` structure.

-----------------------------

.. code-block:: matlab

    B = bg_sumstrexp(t,param,lambda)

Computes the background model ``B`` for a given dipolar pathway amplitude ``lambda``.

