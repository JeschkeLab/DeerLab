.. highlight:: python
.. _dd_randcoil:

***********************
:mod:`dd_randcoil`
***********************

.. autofunction:: deerlab.dd_models.dd_randcoil


Model
=========================================

.. image:: ../images/model_scheme_dd_randcoil.png
   :width: 35%

:math:`P(r) = \frac{3}{(2\pi\nu_0)^{3/2}}4\pi r^2\exp(-\frac{3 r^2}{\nu_0})`

where :math:`\nu_0 = 3/(12\pi r_0 N \nu)^{3/2}`

============== =========== ======== ======== ======== ==================================
 Variable       Symbol     Default   Lower   Upper       Description
============== =========== ======== ======== ======== ==================================
``param(1)``   :math:`N`      50      2        1000    Number of residues
``param(2)``   :math:`R_0`    0.20    0.10     0.40    Segment length
``param(3)``   :math:`\nu`    0.60    0.33     1.00    Scaling exponent
============== =========== ======== ======== ======== ==================================

Example using default parameters:

.. image:: ../images/model_dd_randcoil.png
   :width: 650px

