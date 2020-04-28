Background models
======================

DeerLab includes a collection of parametric models that can be used to model the background signal, i.e. the signal due to randomly distributed spins in the sample that are not part of the spin-labeled protein or object.

All background model functions start with the prefix ``bg_`` and have the calling syntax. The inputs a time axis vector ``t`` (in microseconds), a parameter vector ``param``, and a modulation amplitude ``lambda`` (between 0 and 1) . The length of ``param``, and the meaning of the elements, depends on the particular model. If ``lambda`` is not provided, it is set to one. The output is a background decay vector ``B``, defined over ``t``.

In the following example, the parameter is a single number and represents the spin concentration:

.. code-block:: matlab

        t = linspace(-1,5,301); % us
        c = 100; % uM
        lambda = 0.5;
        B = bg_hom3d(t,c,lambda);

All models are symmetric with respect to ``t=0``, meaning that they depend only on the magnitude of ``t``, ``|t|``.

Background models fall into two categories, physical and phenomenological.


Physical models
------------------------------------

This category involves models that describe particular distributions of spin labels in space. These models depend on physical parameters such as spin concentration, exclusion distances, and dimensionality.


.. rst-class:: func-list

=============================================  ===================================================================  =========================================================
Model function                                 Description                                                          Parameters
=============================================  ===================================================================  =========================================================
:doc:`bg_hom3d<models/bg_hom3d>`                 Homogeneous distribution in 3D                                      Spin concentration
:doc:`bg_hom3dex<models/bg_hom3dex>`             Homogeneous distribution in 3D with hard-shell excluded volume      Spin concentration, exclusion distance
:doc:`bg_homfractal<models/bg_homfractal>`       Homogeneous distribution in fractal dimensions                      Spin concentration, dimensionality (0-6)
=============================================  ===================================================================  =========================================================

In general, it is preferable to use the physical models above instead of the phenomenological models below.


Phenomenological models
------------------------------------

This category contains phenomenological models that represent various mathematical functions that are intended to mimick the background decay, without reference to a particular spatial distribution of spins. The parameters of these models do no have direct physical meaning.

.. rst-class:: func-list

=============================================  ============================================================
Model function                                 Description
=============================================  ============================================================
:doc:`bg_exp<models/bg_exp>`                   Exponential decay
:doc:`bg_strexp<models/bg_strexp>`             Stretched exponential decay
:doc:`bg_sumstrexp<models/bg_sumstrexp>`       Sum of two stretched exponential decays
:doc:`bg_prodstrexp<models/bg_prodstrexp>`     Product of two stretched exponential decays
:doc:`bg_poly1<models/bg_poly1>`               Linear polynomial
:doc:`bg_poly2<models/bg_poly2>`               Quadratic polynomial
:doc:`bg_poly3<models/bg_poly3>`               Cubic polynomial
=============================================  ============================================================

Although the ``bg_exp`` model has the same functional form as ``bg_hom3d``, it is distinct since its parameter is a decay rate constant and not a spin concentration like for ``bg_hom3d``. A similar situation holds for ``bg_strexp`` and ``bg_homfractal``.

.. toctree::
    :maxdepth: 0
    :hidden:
    :includehidden:
    :glob:
    :caption: Background models

    ./models/bg_*