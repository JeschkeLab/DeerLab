.. highlight:: matlab
.. _rd_fourgaussian:


************************
:mod:`rd_fourgaussian`
************************

Sum of four Gaussian distributions parametric model

-----------------------------


Syntax
=========================================

.. code-block:: matlab

        info = rd_fourgaussian()
        P = rd_fourgaussian(r,param)

Parameters
    *   ``t`` - Time axis (N-array)
    *   ``param`` - Model parameters
Returns
    *   ``B`` - Model background (N-array)
    *   ``info`` - Model information (struct)


-----------------------------

Model
=========================================

:math:`P(r) = A_1\sqrt{\frac{2}{\pi}}\frac{1}{\Gamma_1}\exp\left(-\frac{(r-\left<r_1\right>)^2}{\Gamma_1^2}\right) + A_2\sqrt{\frac{2}{\pi}}\frac{1}{\Gamma_2}\exp\left(-\frac{(r-\left<r_2\right>)^2}{\Gamma_2^2}\right) + A_3\sqrt{\frac{2}{\pi}}\frac{1}{\Gamma_3}\exp\left(-\frac{(r-\left<r_3\right>)^2}{\Gamma_3^2}\right) +  (1 - A_1 - A_2 - A_3)\sqrt{\frac{2}{\pi}}\frac{1}{\Gamma_4}\exp\left(-\frac{(r-\left<r_4\right>)^2}{\Gamma_4^2}\right)`

with :math:`\Gamma_i = w_i/\sqrt{2ln(2)}`

========== ======================== ========= ======== ========= ===================================
 Variable   Symbol                    Default   Lower    Upper       Description
========== ======================== ========= ======== ========= ===================================
param(1)   :math:`\left<r_1\right>`     2.5     1.0        20         1st Gaussian mean distance
param(2)   :math:`w_1`                  0.5     0.2        5          1st Gaussian FWHM
param(3)   :math:`\left<r_2\right>`     3.0     1.0        20         2nd Gaussian mean distance
param(4)   :math:`w_2`                  0.5     0.2        5          2nd Gaussian FWHM
param(5)   :math:`\left<r_3\right>`     4.0     1.0        20         3rd Gaussian mean distance
param(6)   :math:`w_3`                  0.5     0.2        5          3rd Gaussian FWHM
param(7)   :math:`\left<r_4\right>`     5.0     1.0        20         4th Gaussian mean distance
param(8)   :math:`w_4`                  0.5     0.2        5          4th Gaussian FWHM
param(9)   :math:`A_1`                  0.25     0          1          1st Gaussian relative amplitude
param(10)  :math:`A_2`                  0.25     0          1          2nd Gaussian relative amplitude
param(10)  :math:`A_3`                  0.25     0          1          3rd Gaussian relative amplitude
========== ======================== ========= ======== ========= ===================================

-----------------------------


Description
=========================================

.. code-block:: matlab

        info = rd_fourgaussian()

Returns an ``info`` structure containing the specifics of the model:

* ``info.model`` -  Full name of the parametric model.
* ``info.nparam`` -  Total number of adjustable parameters.
* ``info.parameters`` - Structure array with information on individual parameters.

-----------------------------


.. code-block:: matlab

    P = rd_fourgaussian(r,param)

Computes the distance distribution model ``P`` from the axis ``r`` according to the parameters array ``param``. The required parameters can also be found in the ``info`` structure.

