Background models
======================

DeerLab includes a collection of parametric models that can be used to model the background signal, i.e. the signal due to randomly distributed spins in the sample that are not part of the spin-labeled protein or object.

All background model functions start with the prefix ``bg_``. They take a time axis vector ``t`` (in microseconds) and a parameter vector ``param`` as inputs and return a background function ``B``, defined over ``t``, as output. In addition, they also take pathway amplitudes ``lambda``,which if not specified are by default set to one. Here is an example:

.. code-block:: matlab

        B = bg_strexp(t,param)
        B = bg_strexp(t,param,lambda)


The first class of models involves physical models that represent homogeneous uniform distributions of spin labels in space.

.. rst-class:: func-list

=============================================  ============================================================
Model function                                 Description
=============================================  ============================================================
:doc:`bg_hom3d<models/bg_hom3d>`                 Homogeneous distribution in 3D
:doc:`bg_homfractal<models/bg_homfractal>`       Homogeneous distribution in 0-6 dimensions
=============================================  ============================================================


The second class contains phenomenological models and lack a strict physical basis.

.. rst-class:: func-list

=============================================  ============================================================
Model function                                 Description
=============================================  ============================================================
:doc:`bg_exp<models/bg_exp>`                   Exponential decay
:doc:`bg_strexp<models/bg_strexp>`             Stretched exponential decay
:doc:`bg_exvol<models/bg_exvol>`               Hard-shell excluded-volume model
:doc:`bg_sumstrexp<models/bg_sumstrexp>`       Sum of two stretched exponential decays
:doc:`bg_prodstrexp<models/bg_prodstrexp>`     Product of two stretched exponential decays
:doc:`bg_poly1<models/bg_poly1>`               Linear polynomial
:doc:`bg_poly2<models/bg_poly2>`               Quadratic polynomial
:doc:`bg_poly3<models/bg_poly3>`               Cubic polynomial
=============================================  ============================================================

.. toctree::
    :maxdepth: 0
    :hidden:
    :includehidden:
    :glob:
    :caption: Background models

    ./models/bg_*