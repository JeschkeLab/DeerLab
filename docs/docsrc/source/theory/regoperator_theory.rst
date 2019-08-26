.. highlight:: matlab
.. _regoperator_theory:

***********************
Regularization operator
***********************

For the inversion of DEER data, the solution can be strongly stabilized if a smoothness criterion is imposed upon :math:`P(r)`. Smoothness should be expected from any probability density distribution obtained from a large number of underlying conformations of the systems containing the spins. In particular, for the major application scenario of spin-labeled biomacromolecules, conformation distribution of the label itself introduces a lower limit to the width of features in :math:`P(r)`. The smoothness criterion can be imposed by employing a discrete differential operator matrix :math:`L_n` of :math:`n^{th}`-order as the regularization matrix [1]_

.. math::	L_nP \longleftrightarrow \frac{\partial^n}{\partial r^n}P(r).

By Fourier analysis of a model distance distribution and its discrete derivatives of different order one can find that a higher order of the differential operator enhances the fast-oscillating components in the :math:`L_nP` term, as has been argued before [2]_. Since the :math:`L_nP` term enters the regularization as a penalty, the larger the order of the regularization matrix :math:`L_n`, the larger the damping of the fast-oscillating components leading to smooth solutions.


.. figure:: ../images/regoperator1.svg
    :align: center

Fourier analysis of the differential operator matrices :math:`L_n` of :math:`n^{th}`-order. With increasing derivative order :math:`n`, the  faster-oscillating components are enhanced, whereas slower-oscillating components are suppressed. The different orders of :math:`L_n` are represented by different colors.


References
=========================================

.. [1] L. Fábregas Ibáñez, G. Jeschke, JMR 300, 2019, 28–40
.. [2] G. Huang, S. Noschese, L. Reichel, Linear Algebra Appl. 502, 2016, 41–57