.. highlight:: python
.. _ex_5pdeer:


***********************
:mod:`ex_5pdeer`
***********************

.. autofunction:: deerlab.ex_models.ex_5pdeer

Model
=========================================

.. image:: ../images/model_scheme_ex_5pdeer.png
   :width: 550px


This experiment model has two modulated pathway and an unmodulated contribution. The kernel is 

.. math::
   K(t,r) =
   [\varLambda_0 + \lambda_1 K_0(t-T_0^{(1)},r) + \lambda_2 K_0(t-T_0^{(2)},r)]
   B(t-T_0^{(1)},\lambda_1) B(t - T_0^{(2)},\lambda_2)

where :math:`T_0^{(1)}=0` and :math:`T_0^{(2)}` are the refocusing times of the two modulated dipolar pathways.


============== ======================== ============= ============ ============ ================================================
 Variable        Symbol                  Start Values     Lower        Upper                Description
============== ======================== ============= ============ ============ ================================================
``param[0]``   :math:`\varLambda_0`          0.4            0            1       Unmodulated pathways, amplitude
``param[1]``   :math:`\lambda_1`             0.4            0            1       1st modulated pathway, amplitude
``param[2]``   :math:`\lambda_2`             0.2            0            1       2nd modulated pathway, amplitude
``param[3]``   :math:`T_0^{(2)}`             5.0            0            20      2nd modulated pathway, refocusing time (μs)
============== ======================== ============= ============ ============ ================================================


Example of a simulated signal using start values:

.. image:: ../images/model_ex_5pdeer.png
   :width: 550px
