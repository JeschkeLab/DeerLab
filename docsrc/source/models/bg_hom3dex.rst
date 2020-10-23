.. highlight:: python
.. _bg_hom3dex:

***********************
:mod:`bg_hom3dex`
***********************

.. autofunction:: deerlab.bg_models.bg_hom3dex

Model
=========================================

.. image:: ../images/model_scheme_bg_hom3dex.png
   :width: 350px

This implements a hard-shell excluded-volume model, with pumped spin concentration ``c`` (first parameter, in μM) and distance of closest approach ``R`` (second parameter, in nm).

The expression for this model is

.. math:: B(t) = \mathrm{exp}\left(-\frac{8\pi^2}{9\sqrt{3}}\alpha(R) \lambda c D |t|\right)`

where :math:`c` is the spin concentration (entered in spins/m\ :sup:`3` into this expression) and :math:`D` is the dipolar constant

.. math:: D = \frac{\mu_0}{4\pi}\frac{(g_\mathrm{e}\mu_\mathrm{B})^2}{\hbar}

The function :math:`\alpha(R)` of the exclusion distance :math:`R` captures the excluded-volume effect. It is a smooth function, but doesn't have an analytical representation. For details, see `Kattnig et al, J.Phys.Chem.B 2013, 117, 16542 <https://pubs.acs.org/doi/abs/10.1021/jp408338q>`_.

============== =============== ============= ============= ============= =================================
 Variable         Symbol        Start Value   Lower bound   Upper bound      Description
============== =============== ============= ============= ============= =================================
``param[0]``    :math:`c`              50         0.01          1000          Spin concentration (μM)
``param[1]``    :math:`R`              1          0.1            20           Exclusion distance (nm)
============== =============== ============= ============= ============= =================================