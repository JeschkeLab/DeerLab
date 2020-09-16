.. highlight:: python
.. _dd_gauss:


***********************
:mod:`dd_gauss`
***********************

.. autofunction:: deerlab.dd_models.dd_gauss

Model
=========================================

:math:`P(r) = \sqrt{\frac{2}{\pi}}\frac{1}{\sigma}\exp\left(-\frac{(r-\left<r\right>)^2}{\sigma^2}\right)`

============== ======================== ============= ============= ============= ===========================
 Variable       Symbol                   Start value   Lower bound   Upper bound      Description
============== ======================== ============= ============= ============= ===========================
``param[0]``   :math:`\left<r\right>`     3.5            1.0              20        Mean (nm)
``param[1]``   :math:`\sigma`             0.2            0.05             2.5       Standard deviation (nm)
============== ======================== ============= ============= ============= ===========================


Example using the start values:

.. image:: ../images/model_dd_gauss.png
   :width: 650px
