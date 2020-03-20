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
   :width: 40%


:math:`P(r) = \sqrt{\frac{2}{\pi}}\frac{1}{\sigma}\exp\left(-\frac{(r-\left<r\right>)^2}{\sqrt(2)\sigma^2}\right)\frac{1}{2}\left(1 + erf\left(\frac{(r-\left<r\right>)}{\sqrt{2}\sigma}\right) \right)`

with :math:`\sigma = w/(2\sqrt{2ln(2)})`

============== ======================== ========= ============= ============= ========================
 Variable       Symbol                    Default   Lower bound   Upper bound      Description
============== ======================== ========= ============= ============= ========================
``param(1)``   :math:`\left<r\right>`     3.5     1.0              20         Mean distance
``param(2)``   :math:`w`                  0.5     0.2              5          FWHM
``param(2)``   :math:`\alpha`             5.0     -15              15         Skewness
============== ======================== ========= ============= ============= ========================


Example using default parameters:

.. image:: ../images/model_dd_skewgauss.png
   :width: 40%


-----------------------------


Description
=========================================

.. code-block:: matlab

        info = dd_skewgauss()

Returns an ``info`` structure containing the specifics of the model:

* ``info.model`` -  Full name of the parametric model.
* ``info.nparam`` -  Total number of adjustable parameters.
* ``info.parameters`` - Structure array with information on individual parameters.

-----------------------------


.. code-block:: matlab

    P = dd_skewgauss(r,param)

Computes the distance distribution model ``P`` from the axis ``r`` according to the parameters array ``param``. The required parameters can also be found in the ``info`` structure.

