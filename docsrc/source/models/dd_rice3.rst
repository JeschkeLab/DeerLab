.. highlight:: python
.. _dd_rice3:


***********************
:mod:`dd_rice3`
***********************

.. autofunction:: deerlab.dd_models.dd_rice3


Model
=========================================

:math:`P(r) = a_1 R(r,\nu_1,\sigma_1) + a_2 R(r,\nu_2,\sigma_2) + a_3 R(r,\nu_3,\sigma_3)`

:math:`R(r,\nu,\sigma) = \frac{\nu^{n/2-1}}{\sigma^2}r^{n/2}\exp\left(-\frac{(r^2+\nu^2)}{2\sigma^2}\right)I_{n/2-1}\left(\frac{r\nu}{\sigma^2} \right)`

where :math:`n=3` and :math:`I_{n/2-1}(x)` is the modified Bessel function of the first kind with order :math:`n/2-1`.
This is a three-dimensional non-central chi distribution, the 3D generalization of the 2D Rice distribution.


============== ======================== ============= ============= ============= =======================================
 Variable         Symbol                 Start Value   Lower bound   Upper bound      Description
============== ======================== ============= ============= ============= =======================================
``param[0]``   :math:`\nu_1`                2.5           1.0           10          1st Rician location (nm)
``param[1]``   :math:`\sigma_1`             0.7           0.1           5           1st Rician spread (nm)
``param[2]``   :math:`a_1`                  0.3             0           1           1st Rician amplitude
``param[3]``   :math:`\nu_2`                4.0           1.0           10          2nd Rician location (nm)
``param[4]``   :math:`\sigma_2`             0.7           0.1           5           2nd Rician spread (nm)
``param[5]``   :math:`a_2`                  0.3           0             1           2nd Rician amplitude
``param[6]``   :math:`\nu_3`                5.0           1.0           10          3rd Rician location (nm)
``param[7]``   :math:`\sigma_3`             0.7           0.1           5           3rd Rician spread (nm)
``param[8]``   :math:`a_3`                  0.3           0             1           3rd Rician amplitude
============== ======================== ============= ============= ============= =======================================


Example using start values:

.. image:: ../images/model_dd_rice3.png
   :width: 650px
