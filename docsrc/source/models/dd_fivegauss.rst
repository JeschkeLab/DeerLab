.. highlight:: matlab
.. _dd_fivegauss:


************************
:mod:`dd_fivegauss`
************************

Sum of five Gaussian distributions

-----------------------------


Syntax
=========================================

.. code-block:: matlab

        info = dd_fivegauss()
        P = dd_fivegauss(r,param)

Parameters
    *   ``r`` - Distance axis (N-array)
    *   ``param`` - Model parameters
Returns
    *   ``P`` - Distance distribution (N-array)
    *   ``info`` - Model information (struct)

-----------------------------

Model
=========================================

:math:`P(r) = A_1\sqrt{\frac{2}{\pi}}\frac{1}{\sigma_1}\exp\left(-\frac{(r-\left<r_1\right>)^2}{\sigma_1^2}\right) + A_2\sqrt{\frac{2}{\pi}}\frac{1}{\sigma_2}\exp\left(-\frac{(r-\left<r_2\right>)^2}{\sigma_2^2}\right) + A_3\sqrt{\frac{2}{\pi}}\frac{1}{\sigma_3}\exp\left(-\frac{(r-\left<r_3\right>)^2}{\sigma_3^2}\right) +  A_4\sqrt{\frac{2}{\pi}}\frac{1}{\sigma_4}\exp\left(-\frac{(r-\left<r_4\right>)^2}{\sigma_4^2}\right) + (1 - A_1 - A_2 - A_3 - A_4)\sqrt{\frac{2}{\pi}}\frac{1}{\sigma_5}\exp\left(-\frac{(r-\left<r_5\right>)^2}{\sigma_5^2}\right)`

with :math:`\sigma_i = w_i/\sqrt{2ln(2)}`

============== ======================== ========= ======== ========= ===================================
 Variable       Symbol                    Default   Lower    Upper       Description
============== ======================== ========= ======== ========= ===================================
``param(1)``   :math:`\left<r_1\right>`     2.5     1.0        20         1st Gaussian mean distance
``param(2)``   :math:`w_1`                  0.5     0.2        5          1st Gaussian FWHM
``param(3)``   :math:`\left<r_2\right>`     3.0     1.0        20         2nd Gaussian mean distance
``param(4)``   :math:`w_2`                  0.5     0.2        5          2nd Gaussian FWHM
``param(5)``   :math:`\left<r_3\right>`     3.5     1.0        20         3rd Gaussian mean distance
``param(6)``   :math:`w_3`                  0.5     0.2        5          3rd Gaussian FWHM
``param(7)``   :math:`\left<r_4\right>`     4.5     1.0        20         4th Gaussian mean distance
``param(8)``   :math:`w_4`                  0.5     0.2        5          4th Gaussian FWHM
``param(9)``   :math:`\left<r_5\right>`     5.0     1.0        20         5th Gaussian mean distance
``param(10)``   :math:`w_5`                 0.5     0.2        5          5th Gaussian FWHM
``param(11)``   :math:`A_1`                 0.2     0          1          1st Gaussian relative amplitude
``param(12)``  :math:`A_2`                  0.2     0          1          2nd Gaussian relative amplitude
``param(13)``  :math:`A_3`                  0.2     0          1          3rd Gaussian relative amplitude
``param(14)``  :math:`A_4`                  0.2     0          1          4th Gaussian relative amplitude
============== ======================== ========= ======== ========= ===================================


Example using default parameters:

.. image:: ../images/model_dd_fivegauss.png
   :width: 40%


-----------------------------


Description
=========================================

.. code-block:: matlab

        info = dd_fivegauss()

Returns an ``info`` structure containing the specifics of the model:

* ``info.model`` -  Full name of the parametric model.
* ``info.nparam`` -  Total number of adjustable parameters.
* ``info.parameters`` - Structure array with information on individual parameters.

-----------------------------


.. code-block:: matlab

    P = dd_fivegauss(r,param)

Computes the distance distribution model ``P`` from the axis ``r`` according to the parameters array ``param``. The required parameters can also be found in the ``info`` structure.

