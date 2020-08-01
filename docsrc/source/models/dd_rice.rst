.. highlight:: python
.. _dd_rice:


***********************
:mod:`dd_rice`
***********************


.. autofunction:: deerlab.dd_models.dd_rice


Model
=========================================

:math:`P(r) = \frac{\nu^{n/2-1}}{\sigma^2}r^{n/2}\exp\left(-\frac{(r^2+\nu^2)}{2\sigma^2}\right)I_{n/2-1}\left(\frac{r\nu}{\sigma^2} \right)`

where :math:`n=3` and :math:`I_{n/2-1}(x)` is the modified Bessel function of the first kind with order :math:`n/2-1`. This is a three-dimensional non-central chi distribution, the 3D generalization of the 2D Rice distribution.

============== ======================== ========= ============= ============= ========================
 Variable       Symbol                    Default   Lower bound   Upper bound      Description
============== ======================== ========= ============= ============= ========================
``param(1)``   :math:`\nu`                3.5     1.0              10         Center
``param(2)``   :math:`\sigma`             0.7     0.1              5          Width
============== ======================== ========= ============= ============= ========================


Example using default parameters:

.. image:: ../images/model_dd_rice.png
   :width: 650px

