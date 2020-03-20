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
   :width: 40%

:math:`P(r) = \frac{\beta}{2\sigma\Gamma(1/\beta)}\exp\left(-\left(\frac{(r-\left<r\right>)}{\sigma}\right)^\beta \right)`

with :math:`\sigma = w/(2\sqrt{2ln(2)})`

============== ======================== ========= ============= ============= ========================
 Variable       Symbol                    Default   Lower bound   Upper bound      Description
============== ======================== ========= ============= ============= ========================
``param(1)``   :math:`\left<r\right>`     3.5     1.0              20         Mean distance
``param(2)``   :math:`w`                  0.5     0.2              5          FWHM
``param(2)``   :math:`\beta`              5.0     0.25             15         Kurtosis
============== ======================== ========= ============= ============= ========================


Example using default parameters:

.. image:: ../images/model_dd_gengauss.png
   :width: 40%


-----------------------------


Description
=========================================

.. code-block:: matlab

        info = dd_gengauss()

Returns an ``info`` structure containing the specifics of the model:

* ``info.model`` -  Full name of the parametric model.
* ``info.nparam`` -  Total number of adjustable parameters.
* ``info.parameters`` - Structure array with information on individual parameters.

-----------------------------


.. code-block:: matlab

    P = dd_gengauss(r,param)

Computes the distance distribution model ``P`` from the axis ``r`` according to the parameters array ``param``. The required parameters can also be found in the ``info`` structure.

