Background models
======================

DeerLab includes a collection of parametric models that can be used to model the background signal, i.e. the signal due to randomly distributed spins in the sample that are not part of the spin-labeled protein or object.

All background model functions start with the prefix ``bg_``. They take a time axis vector ``t`` (in microseconds) and a parameter vector ``param`` as inputs and return a background function ``B``, defined over ``t``, as output. Here is an example:

.. code-block:: matlab

        B = bg_strexp(t,param)



The first class of models involves exponential decays of various forms. The exponential and stretched exponential decays are physical models that represent homogeneous uniform distributions of spin labels.

.. rst-class:: func-list

=============================================  ============================================================
Model function                                 Description
=============================================  ============================================================
:doc:`bg_exp<models/bg_exp>`                   Exponential decay
:doc:`bg_strexp<models/bg_strexp>`             Stretched exponential decay
:doc:`bg_sumstrexp<models/bg_sumstrexp>`       Sum of two stretched exponential decays
:doc:`bg_prodstrexp<models/bg_prodstrexp>`     Product of two stretched exponential decays
=============================================  ============================================================

The second class contains polynomial functions. These are purely phenomenological and lack a physical basis.

.. rst-class:: func-list

=============================================  ============================================================
Model function                                 Description
=============================================  ============================================================
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