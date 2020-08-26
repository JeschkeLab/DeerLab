.. highlight:: python
.. _dd_wormchain:

***********************
:mod:`dd_wormchain`
***********************

.. autofunction:: deerlab.dd_models.dd_wormchain

Model
=========================================


============== ============ ============= ============= ============= =========================
 Variable         Symbol     Start Value   Lower bound   Upper bound      Description
============== ============ ============= ============= ============= =========================
``param[0]``   :math:`L`      3.7            1.5            10         Contour length (nm)
``param[1]``   :math:`L_p`    10             2              100        Persistence length (nm)
============== ============ ============= ============= ============= =========================

Example using start values:

.. image:: ../images/model_dd_wormchain.png
   :width: 650px


References
=========================================

[1] J. Wilhelm, E. Frey, Phys. Rev. Lett. 77(12), 2581-2584 (1996)
DOI:  `10.1103/PhysRevLett.77.2581 <https://doi.org/10.1103/PhysRevLett.77.2581>`_