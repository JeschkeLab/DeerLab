.. highlight:: python
.. _dd_skewgauss:


***********************
:mod:`dd_skewgauss`
***********************

.. autofunction:: deerlab.dd_models.dd_skewgauss

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
