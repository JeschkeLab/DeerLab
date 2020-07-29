.. highlight:: matlab
.. _dd_cos:


***********************
:mod:`dd_cos`
***********************

Raised-cosine distribution

-----------------------------


Syntax
=========================================

.. code-block:: matlab

        info = dd_cos()
        P = dd_cos(r,param)

Parameters
    *   ``r`` - Distance axis (N-array), in nanometers
    *   ``param`` - Model parameters
Returns
    *   ``P`` - Distance distribution (N-array)
    *   ``info`` - Model information (struct)

-----------------------------

Model
=========================================


This provides a `raised-cosine distribution <https://en.wikipedia.org/wiki/Raised_cosine_distribution>`_, defined by 
:math:`P(r) = \frac{1}{2w}\cos\left(\frac{r-r_0}{w}\pi\right)` for :math:`r_0-w \le r \le r_0+w`, and zero otherwise.

============== ======================== ========= ============= ============= ========================
 Variable       Symbol                    Default   Lower bound   Upper bound      Description
============== ======================== ========= ============= ============= ========================
``param(1)``   :math:`r_0`                 3.0       0.1              20          center, in nm
``param(2)``   :math:`w`                   0.5       0.1               5          fwhm, in nm
============== ======================== ========= ============= ============= ========================


Example using default parameters:

.. image:: ../images/model_dd_cos.png
   :width: 650px


-----------------------------


Description
=========================================

.. code-block:: matlab

        info = dd_cos()

Returns an ``info`` structure containing the information of the model parameters and boundaries.

* ``info(n).Index`` -  Index of the parameter in the ``param`` array.
* ``info(n).Parameter`` -  Description of the n-th parameter.
* ``info(n).Lower`` -  Lower bound of the n-th parameter.
* ``info(n).Upper`` -  Upper bound of the n-th parameter.
* ``info(n).Start`` -  Start value of the n-th parameter.

-----------------------------


.. code-block:: matlab

    P = dd_cos(r,[3 0.5]])

Computes the distance distribution model ``P`` from the axis ``r`` according to the parameters array ``param``. The required parameters can also be found in the ``info`` structure.

