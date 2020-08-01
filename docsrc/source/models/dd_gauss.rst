.. highlight:: python
.. _dd_gauss:


***********************
:mod:`dd_gauss`
***********************

.. autofunction:: deerlab.dd_models.dd_gauss

Model
=========================================

:math:`P(r) = \sqrt{\frac{2}{\pi}}\frac{1}{\sigma}\exp\left(-\frac{(r-\left<r\right>)^2}{\sigma^2}\right)`

with :math:`\sigma = \mathrm{FWHM}/\sqrt{2ln(2)}`

============== ======================== ========= ============= ============= ========================
 Variable       Symbol                    Default   Lower bound   Upper bound      Description
============== ======================== ========= ============= ============= ========================
``param(1)``   :math:`r0`                 3.5         1.0              20         center distance (nm)
``param(2)``   :math:`\mathrm{FWHM}`      0.5         0.2              5          FWHM (nm)
============== ======================== ========= ============= ============= ========================


Example using default parameters:

.. image:: ../images/model_dd_gauss.png
   :width: 650px
