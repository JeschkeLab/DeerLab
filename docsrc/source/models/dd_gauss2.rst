.. highlight:: matlab
.. _dd_gauss2:


***********************
:mod:`dd_gauss2`
***********************

Sum of two Gaussian distributions

-----------------------------


Syntax
=========================================

.. code-block:: matlab

        info = dd_gauss2()
        P = dd_gauss2(r,param)

Parameters
    *   ``r`` - Distance axis (N-array)
    *   ``param`` - Model parameters
Returns
    *   ``P`` - Distance distribution (N-array)
    *   ``info`` - Model information (struct)


-----------------------------

Model
=========================================

:math:`P(r) = a_1\sqrt{\frac{2}{\pi}}\frac{1}{\sigma_1}\exp\left(-\frac{(r-\left<r_1\right>)^2}{\sigma_1^2}\right) + a_2\sqrt{\frac{2}{\pi}}\frac{1}{\sigma_2}\exp\left(-\frac{(r-\left<r_2\right>)^2}{\sigma_2^2}\right)`

with :math:`\sigma_i = \mathrm{FWHM}_i/\sqrt{2ln(2)}`


============== ======================== ========= ======== ========= ===================================
 Variable       Symbol                    Default   Lower    Upper       Description
============== ======================== ========= ======== ========= ===================================
``param(1)``   :math:`\left<r_1\right>`     2.5     1.0        20       1st Gaussian center distance
``param(2)``   :math:`\mathrm{FWHM}_1`      0.5     0.2        5        1st Gaussian FWHM
``param(3)``   :math:`a_1`                  0.5     0          1        1st Gaussian amplitude
``param(4)``   :math:`\left<r_2\right>`     3.5     1.0        20       2nd Gaussian center distance
``param(5)``   :math:`\mathrm{FWHM}_2`      0.5     0.2        5        2nd Gaussian FWHM
``param(6)``   :math:`a_2`                  0.5     0          1        2nd Gaussian amplitude
============== ======================== ========= ======== ========= ===================================

Example using default parameters:

.. image:: ../images/model_dd_gauss2.png
   :width: 650px


-----------------------------


Description
=========================================

.. code-block:: matlab

        info = dd_gauss2()

Returns an ``info`` structure containing the information of the model parameters and boundaries.

* ``info(n).Index`` -  Index of the parameter in the ``param`` array.
* ``info(n).Parameter`` -  Description of the n-th parameter.
* ``info(n).Lower`` -  Lower bound of the n-th parameter.
* ``info(n).Upper`` -  Upper bound of the n-th parameter.
* ``info(n).Start`` -  Start value of the n-th parameter.

-----------------------------


.. code-block:: matlab

    P = dd_gauss2(r,param)

Computes the distance distribution model ``P`` from the axis ``r`` according to the parameters array ``param``. The required parameters can also be found in the ``info`` structure.

