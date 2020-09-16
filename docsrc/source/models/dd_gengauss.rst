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

============== ========================== ============= ============= ============= =======================================
 Variable         Symbol                   Start Value   Lower bound   Upper bound      Description
============== ========================== ============= ============= ============= =======================================
``param[0]``   :math:`r_0`                     3.5          1.0              20         Mean (nm)
``param[1]``   :math:`\sigma`                  0.2          0.05             2.5        Spread (nm)
``param[2]``   :math:`\beta`                   5.0          0.25             15         kurtosis
============== ========================== ============= ============= ============= =======================================


Example using start values:

.. image:: ../images/model_dd_gengauss.png
   :width: 650px
