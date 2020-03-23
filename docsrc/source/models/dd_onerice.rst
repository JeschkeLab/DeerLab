.. highlight:: matlab
.. _dd_onerice:


***********************
:mod:`dd_onerice`
***********************

Rician distribution parametric model

-----------------------------


Syntax
=========================================

.. code-block:: matlab

        info = dd_onerice()
        P = dd_onerice(r,param)

Parameters
    *   ``r`` - Distance axis (N-array)
    *   ``param`` - Model parameters
Returns
    *   ``P`` - Distance distribution (N-array)
    *   ``info`` - Model information (struct)

-----------------------------

Model
=========================================

:math:`P(r) = \frac{\left<r\right>^{L-1}}{\sigma^2}r^L\exp\left(-\frac{(r^2+\left<r\right>^2)}{2\sigma^2}\right)I_{L-1}\left(\frac{r\left<r\right>}{\sigma^2} \right)`

where :math:`I_0(x)` is the modified Bessel function of the first kind with order zero.

============== ======================== ========= ============= ============= ========================
 Variable       Symbol                    Default   Lower bound   Upper bound      Description
============== ======================== ========= ============= ============= ========================
``param(1)``   :math:`\left<r\right>`     3.5     1.0              10         Mean distance
``param(2)``   :math:`\sigma`             0.7     0.1              5          Standard deviation
============== ======================== ========= ============= ============= ========================


Example using default parameters:

.. image:: ../images/model_dd_onerice.png
   :width: 40%


-----------------------------


Description
=========================================

.. code-block:: matlab

        info = dd_onerice()

Returns an ``info`` structure containing the specifics of the model:

* ``info.model`` -  Full name of the parametric model.
* ``info.nparam`` -  Total number of adjustable parameters.
* ``info.parameters`` - Structure array with information on individual parameters.

-----------------------------


.. code-block:: matlab

    P = dd_onerice(r,param)

Computes the distance distribution model ``P`` from the axis ``r`` according to the parameters array ``param``. The required parameters can also be found in the ``info`` structure.

