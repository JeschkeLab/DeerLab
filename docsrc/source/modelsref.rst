.. _modelsref:

Models
======================

.. _modelsref_dd:


Distribution models
------------------------------------

Basis functions
........................................
This group of distribution models is based on function that can serve as basis functions for modeling general distributions.

.. rst-class:: func-list

============================= ===================================================
  Model function                 Description
============================= ===================================================
:ref:`dd_rice`                  Single 3D-Rice distribution
:ref:`dd_rice2`                 Two 3D-Rice distributions
:ref:`dd_rice3`                 Three 3D-Rice distributions
:ref:`dd_rice4`                 Four 3D-Rice distributions
:ref:`dd_rice5`                 Five 3D-Rice distributions
:ref:`dd_gauss`                 Single Gaussian distribution
:ref:`dd_gauss2`                Two Gaussians distributions
:ref:`dd_gauss3`                Three Gaussians distributions
:ref:`dd_gauss4`                Four Gaussians distributions
:ref:`dd_gauss5`                Five Gaussians distributions
:ref:`dd_gengauss`              Single generalized Gaussian distribution
:ref:`dd_skewgauss`             Single skew Gaussian distribution
:ref:`dd_cos`                   Raised-cosine distribution
============================= ===================================================

Disordered physical models
........................................

This group of distribution models represents spin labels attached to disordered macromolecular structures.

.. rst-class:: func-list

=========================== ===================================================
  Model function                                       Description
=========================== ===================================================
:ref:`dd_wormchain`           Worm-like chain
:ref:`dd_wormgauss`           Worm-like chain with Gaussian convolution
:ref:`dd_randcoil`            Random coil
=========================== ===================================================


Sphere/shell models
........................................

This group of distribution models represents spin labels in simple partitions of 3D space.

.. rst-class:: func-list

============================== =============================================================================
  Model function                                       Description
============================== =============================================================================
:ref:`dd_sphere`                 Labels distributed in a sphere
:ref:`dd_spheresurf`             Labels distributed on a sphere's surface
:ref:`dd_spherepoint`            One label fixed, other label distributed in a sphere
:ref:`dd_shell`                  Labels distributed on a spherical shell
:ref:`dd_shellsphere`            Labels distributed on a sphere inside a spherical shell
:ref:`dd_shellshell`             Labels distributed on a spherical shell inside another spherical shell
:ref:`dd_shellvoidshell`         Labels distributed on two spherical shells separated by a void
:ref:`dd_shellvoidsphere`        Labels distributed on a sphere inside a spherical shell separated by a void 
============================== =============================================================================


Toy models
........................................

This group contains distribution models that have absolutely no physical relevance. They are useful as distributions to test numerical algorithms.

.. rst-class:: func-list

============================== =============================================================================
  Model function                                       Description
============================== =============================================================================
:ref:`dd_circle`                Semi-circle distribution
:ref:`dd_triangle`              Triangle distribution
:ref:`dd_uniform`               Uniform distribution
============================== =============================================================================



.. _modelsref_bg:

Background models
------------------------------------

Physical models
...........................................

This category involves models that describe particular distributions of spin labels in space. These models depend on physical parameters such as spin concentration, exclusion distances, and dimensionality.

.. rst-class:: func-list

=============================================  =================================================================== 
Model function                                 Description                                                         
=============================================  ===================================================================
:doc:`bg_hom3d<models/bg_hom3d>`                 Homogeneous distribution in 3D
:doc:`bg_hom3dex<models/bg_hom3dex>`             Homogeneous distribution in 3D with hard-shell excluded volume
:doc:`bg_homfractal<models/bg_homfractal>`       Homogeneous distribution in fractal dimensions
=============================================  ===================================================================

Phenomenological models
...........................................

This category comprises phenomenological models that represent various mathematical functions that are intended to mimick the background decay, without reference to a particular spatial distribution of spins. The parameters of these models do no have direct physical meaning.

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


.. _modelsref_ex:


Experiment models
------------------------------------

.. rst-class:: func-list

=========================================== ===============================================================
Model function                              Description
=========================================== ===============================================================
:doc:`ex_4pdeer<models/ex_4pdeer>`          4-pulse DEER (single modulated pathway)
:doc:`ex_ovl4pdeer<models/ex_ovl4pdeer>`     4-pulse DEER (including 2+1 pathway) 
:doc:`ex_5pdeer<models/ex_5pdeer>`          5-pulse DEER
:doc:`ex_7pdeer<models/ex_7pdeer>`           7-pulse DEER
=========================================== ===============================================================



.. toctree::
    :maxdepth: 0
    :hidden:
    :includehidden:
    :glob:
    :caption: Models

    ./models/dd_*
    ./models/bg_*
    ./models/ex_*