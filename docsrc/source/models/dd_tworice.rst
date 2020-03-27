.. highlight:: matlab
.. _dd_tworice:

***********************
:mod:`dd_tworice`
***********************

Sum of two rician distributions parametric model

-----------------------------


Syntax
=========================================

.. code-block:: matlab

        info = dd_tworice()
        P = dd_tworice(r,param)

Parameters
    *   ``r`` - Distance axis (N-array)
    *   ``param`` - Model parameters
Returns
    *   ``P`` - Distance distribution (N-array)
    *   ``info`` - Model information (struct)


-----------------------------

Model
=========================================

:math:`P(r) = A_1 R(r,\nu_1,\sigma_1) + (1-A_1) R(r,\nu_2,\sigma_2)`

:math:`R(r,\nu,\sigma) = \frac{\nu^{n/2-1}}{\sigma^2}r^{n/2}\exp\left(-\frac{(r^2+\nu^2)}{2\sigma^2}\right)I_{n/2-1}\left(\frac{r\nu}{\sigma^2} \right)`

where :math:`n=3` and :math:`I_{n/2-1}(x)` is the modified Bessel function of the first kind with order :math:`n/2-1`.
This is a three-dimensional non-central chi distribution, the 3D generalization of the 2D Rice distribution.

============== ======================== ========= ======== ======== ===============================
 Variable       Symbol                    Default   Lower   Upper       Description
============== ======================== ========= ======== ======== ===============================
``param(1)``   :math:`\nu_1`                2.5     1.0      10      1st mean distance
``param(2)``   :math:`\sigma_1`             0.4     0.1      5       1st standard deviation
``param(3)``   :math:`\nu_2`                4.0     1.0      10      2nd mean distance
``param(4)``   :math:`\sigma_2`             0.4     0.1      5       2nd standard deviation
``param(5)``   :math:`A_1`                  0.5     0        1       fraction of 1st component
============== ======================== ========= ======== ======== ===============================


Example using default parameters:

.. image:: ../images/model_dd_tworice.png
   :width: 40%


-----------------------------


Description
=========================================

.. code-block:: matlab

        info = dd_tworice()

Returns an ``info`` structure containing the specifics of the model:

* ``info.model`` -  Full name of the parametric model.
* ``info.nparam`` -  Total number of adjustable parameters.
* ``info.parameters`` - Structure array with information on individual parameters.

-----------------------------


.. code-block:: matlab

    P = dd_tworice(r,param)

Computes the distance distribution model ``P`` from the axis ``r`` according to the parameters array ``param``. The required parameters can also be found in the ``info`` structure.

