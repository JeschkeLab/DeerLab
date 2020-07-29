.. highlight:: matlab
.. _dd_gengauss:


***********************
:mod:`dd_gengauss`
***********************

Generalized Gaussian distribution parametric model

-----------------------------


Syntax
=========================================

.. code-block:: matlab

        info = dd_gengauss()
        P = dd_gengauss(r,param)

Parameters
    *   ``r`` - Distance axis (N-array)
    *   ``param`` - Model parameters
Returns
    *   ``P`` - Distance distribution (N-array)
    *   ``info`` - Model information (struct)

-----------------------------

Model
=========================================

.. image:: ../images/model_scheme_dd_gengauss.png
   :width: 650px

:math:`P(r) = \frac{\beta}{2\sigma\Gamma(1/\beta)}\exp\left(-\left(\frac{(r-\left<r\right>)}{\sigma}\right)^\beta \right)`

with :math:`\sigma = w/(2\sqrt{2ln(2)})`

============== ======================== ========= ============= ============= ========================
 Variable       Symbol                    Default   Lower bound   Upper bound      Description
============== ======================== ========= ============= ============= ========================
``param(1)``   :math:`r_0`                3.5     1.0              20         center distance (nm)
``param(2)``   :math:`w`                  0.5     0.2              5          FWHM (nm)
``param(2)``   :math:`\beta`              5.0     0.25             15         kurtosis
============== ======================== ========= ============= ============= ========================


Example using default parameters:

.. image:: ../images/model_dd_gengauss.png
   :width: 650px


-----------------------------


Description
=========================================

.. code-block:: matlab

        info = dd_gengauss()

Returns an ``info`` structure containing the information of the model parameters and boundaries.

* ``info(n).Index`` -  Index of the parameter in the ``param`` array.
* ``info(n).Parameter`` -  Description of the n-th parameter.
* ``info(n).Lower`` -  Lower bound of the n-th parameter.
* ``info(n).Upper`` -  Upper bound of the n-th parameter.
* ``info(n).Start`` -  Start value of the n-th parameter.

-----------------------------


.. code-block:: matlab

    P = dd_gengauss(r,param)

Computes the distance distribution model ``P`` from the axis ``r`` according to the parameters array ``param``. The required parameters can also be found in the ``info`` structure.

