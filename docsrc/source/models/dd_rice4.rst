.. highlight:: matlab
.. _dd_rice4:


***********************
:mod:`dd_rice4`
***********************

Sum of four 3D-Rice distributions

-----------------------------


Syntax
=========================================

.. code-block:: matlab

        info = dd_rice4()
        P = dd_rice4(r,param)

Parameters
    *   ``r`` - Distance axis (N-array)
    *   ``param`` - Model parameters
Returns
    *   ``P`` - Distance distribution (N-array)
    *   ``info`` - Model information (struct)


-----------------------------

Model
=========================================

:math:`P(r) = a_1 R(r,\nu_1,\sigma_1) + a_2 R(r,\nu_2,\sigma_2) + a_3 R(r,\nu_3,\sigma_3) + a_4 R(r,\nu_4,\sigma_4)`

:math:`R(r,\nu,\sigma) = \frac{\nu^{n/2-1}}{\sigma^2}r^{n/2}\exp\left(-\frac{(r^2+\nu^2)}{2\sigma^2}\right)I_{n/2-1}\left(\frac{r\nu}{\sigma^2} \right)`

where :math:`n=4` and :math:`I_{n/2-1}(x)` is the modified Bessel function of the first kind with order :math:`n/2-1`.
This is a three-dimensional non-central chi distribution, the 3D generalization of the 2D Rice distribution.

============== ======================== ========= ======== ======== ===============================
 Variable       Symbol                    Default   Lower   Upper       Description
============== ======================== ========= ======== ======== ===============================
``param(1)``   :math:`\nu_1`                2.5     1.0      10      1st Rician center distance
``param(2)``   :math:`\sigma_1`             0.7     0.1      5       1st Rician width
``param(3)``   :math:`a_1`                  0.25     0       1       1st Rician amplitude
``param(4)``   :math:`\nu_2`                3.5     1.0      10      2nd Rician center distance
``param(5)``   :math:`\sigma_2`             0.7     0.1      5       2nd Rician width
``param(6)``   :math:`a_2`                  0.25     0       1       2nd Rician amplitude
``param(7)``   :math:`\nu_3`                4.5     1.0      10      3rd Rician center distance
``param(8)``   :math:`\sigma_3`             0.7     0.1      5       3rd Rician width
``param(9)``   :math:`a_3`                  0.25     0       1       3rd Rician amplitude
``param(10)``  :math:`\nu_4`                5.5     1.0      10      4th Rician center distance
``param(11)``  :math:`\sigma_4`             0.7     0.1      5       4th Rician width
``param(12)``  :math:`a_4`                  0.25     0       1       4th Rician amplitude
============== ======================== ========= ======== ======== ===============================


Example using default parameters:

.. image:: ../images/model_dd_rice4.png
   :width: 550px

-----------------------------


Description
=========================================

.. code-block:: matlab

        info = dd_rice4()

Returns an ``info`` structure containing the information of the model parameters and boundaries.

* ``info(n).Index`` -  Index of the parameter in the ``param`` array.
* ``info(n).Parameter`` -  Description of the n-th parameter.
* ``info(n).Lower`` -  Lower bound of the n-th parameter.
* ``info(n).Upper`` -  Upper bound of the n-th parameter.
* ``info(n).Start`` -  Start value of the n-th parameter.

-----------------------------


.. code-block:: matlab

    P = dd_rice4(r,param)

Computes the distance distribution model ``P`` from the axis ``r`` according to the parameters array ``param``. The required parameters can also be found in the ``info`` structure.

