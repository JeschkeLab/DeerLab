.. highlight:: matlab
.. _dd_skewgauss:


***********************
:mod:`dd_skewgauss`
***********************

Skew Gaussian distribution parametric model

-----------------------------


Syntax
=========================================

.. code-block:: matlab

        info = dd_skewgauss()
        P = dd_skewgauss(r,param)

Parameters
    *   ``r`` - Distance axis (N-array)
    *   ``param`` - Model parameters
Returns
    *   ``P`` - Distance distribution (N-array)
    *   ``info`` - Model information (struct)

-----------------------------

Model
=========================================

.. image:: ../images/model_scheme_dd_skewgauss.png
   :width: 650px


:math:`P(r) = \sqrt{\frac{2}{\pi}}\frac{1}{\sigma}\exp\left(-\frac{(r-\left<r\right>)^2}{\sqrt(2)\sigma^2}\right)\frac{1}{2}\left(1 + erf\left(\frac{(r-\left<r\right>)}{\sqrt{2}\sigma}\right) \right)`

with :math:`\sigma = w/(2\sqrt{2ln(2)})`

============== ======================== ========= ============= ============= ========================
 Variable       Symbol                    Default   Lower bound   Upper bound      Description
============== ======================== ========= ============= ============= ========================
``param(1)``   :math:`r_0`                3.5     1.0              20         Center distance (nm)
``param(2)``   :math:`w`                  0.5     0.2              5          FWHM (nm)
``param(2)``   :math:`\alpha`             5.0     -15              15         Skewness
============== ======================== ========= ============= ============= ========================


Example using default parameters:

.. image:: ../images/model_dd_skewgauss.png
   :width: 650px


-----------------------------


Description
=========================================

.. code-block:: matlab

        info = dd_skewgauss()

Returns an ``info`` structure containing the information of the model parameters and boundaries.

* ``info(n).Index`` -  Index of the parameter in the ``param`` array.
* ``info(n).Parameter`` -  Description of the n-th parameter.
* ``info(n).Lower`` -  Lower bound of the n-th parameter.
* ``info(n).Upper`` -  Upper bound of the n-th parameter.
* ``info(n).Start`` -  Start value of the n-th parameter.

-----------------------------


.. code-block:: matlab

    P = dd_skewgauss(r,param)

Computes the distance distribution model ``P`` from the axis ``r`` according to the parameters array ``param``. The required parameters can also be found in the ``info`` structure.

