.. highlight:: matlab
.. _td_prodstrexp:


***********************
:mod:`td_prodstrexp`
***********************

Product of two stretched exponentials background parametric model

-----------------------------


Syntax
=========================================

.. code-block:: matlab

        info = td_prodstrexp()
        P = td_prodstrexp(r,param)

Parameters
    *   ``t`` - Time axis (N-array)
    *   ``param`` - Model parameters
Returns
    *   ``B`` - Model background (N-array)
    *   ``info`` - Model information (struct)



-----------------------------

Model
=========================================

:math:`B(t) = e^{-(k_1t)^{d_1/3}}e^{-(k_2t)^{d_2/3}}`

========== ============= ========= ============= ============= ==============================
 Variable   Symbol        Default   Lower bound   Upper bound      Description
========== ============= ========= ============= ============= ==============================
param(1)   :math:`k_1`      3.5         0            200         1st strexp decay rate
param(2)   :math:`d_1`      3           0            6           1st strexp fractal dimension
param(3)   :math:`k_2`      3.5         0            200         2nd strexp decay rate
param(4)   :math:`d_2`      3           0            6           2nd strexp fractal dimension
========== ============= ========= ============= ============= ==============================

-----------------------------


Description
=========================================

.. code-block:: matlab

        info = td_prodstrexp()

Returns an ``info`` structure containing the specifics of the model:

* ``info.model`` -  Full name of the parametric model.
* ``info.nparam`` -  Total number of adjustable parameters.
* ``info.parameters`` - Structure array with information on individual parameters.

-----------------------------


.. code-block:: matlab

    B = td_prodstrexp(t,param)

Computes the background model ``B`` from the axis ``t`` according to the parameters array ``param``. The required parameters can also be found in the ``info`` structure.

