.. highlight:: matlab
.. _td_strexp:


***********************
:mod:`td_strexp`
***********************

Stretched exponential background parametric model

-----------------------------


Syntax
=========================================

.. code-block:: matlab

        info = td_strexp()
        P = td_strexp(r,param)

Parameters
    *   ``t`` - Time axis (N-array)
    *   ``param`` - Model parameters
Returns
    *   ``B`` - Model background (N-array)
    *   ``info`` - Model information (struct)


-----------------------------

Model
=========================================

:math:`B(t) = e^{-(kt)^{d/3}}`

========== ========== ========= ============= ============= ========================
 Variable   Symbol     Default   Lower bound   Upper bound      Description
========== ========== ========= ============= ============= ========================
param(1)   :math:`k`      3.5      0              200           Decay rate
param(2)   :math:`d`      3        0              6             Fractal dimension
========== ========== ========= ============= ============= ========================

-----------------------------


Description
=========================================

.. code-block:: matlab

        info = td_strexp()

Returns an ``info`` structure containing the specifics of the model:

* ``info.model`` -  Full name of the parametric model.
* ``info.nparam`` -  Total number of adjustable parameters.
* ``info.parameters`` - Structure array with information on individual parameters.

-----------------------------


.. code-block:: matlab

    B = td_strexp(t,param)

Computes the background model ``B`` from the axis ``t`` according to the parameters array ``param``. The required parameters can also be found in the ``info`` structure.

