.. highlight:: python
.. _dd_gengauss:


***********************
:mod:`dd_gengauss`
***********************


.. autofunction:: deerlab.dd_models.dd_gengauss

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
