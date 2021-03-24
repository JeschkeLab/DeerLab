.. _modelsref:

Models
======================

.. _modelsref_dd:


Distribution models
------------------------------------

Basis functions
........................................

This group of distribution models is based on functions that can serve as basis functions for modeling more general distributions.

.. rst-class:: func-list

============================== =============================================================================
  Model function                 Description
============================== =============================================================================
:ref:`dd_rice`                  3D-Rice distribution
:ref:`dd_rice2`                 Two 3D-Rice distributions
:ref:`dd_rice3`                 Three 3D-Rice distributions
:ref:`dd_gauss`                 Gaussian distribution
:ref:`dd_gauss2`                Two Gaussian distributions
:ref:`dd_gauss3`                Three Gaussian distributions
:ref:`dd_gengauss`              Generalized Gaussian distribution
:ref:`dd_skewgauss`             Skew Gaussian distribution
:ref:`dd_cos`                   Raised-cosine distribution
============================== =============================================================================

Disordered physical models
........................................

This group of distribution models represents spin labels attached to disordered macromolecular structures.

.. rst-class:: func-list

============================== =============================================================================
  Model function                                       Description
============================== =============================================================================
:ref:`dd_wormchain`              Worm-like chain
:ref:`dd_wormgauss`              Worm-like chain with Gaussian convolution
:ref:`dd_randcoil`               Random coil
============================== =============================================================================


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

This group contains distribution models that have absolutely no physical relevance. They are useful for testing numerical algorithms.

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

==============================  ============================================================================= 
Model function                                Description                                                         
==============================  =============================================================================
:ref:`bg_hom3d`                   Homogeneous distribution in 3D
:ref:`bg_hom3dex`                 Homogeneous distribution in 3D with hard-shell excluded volume
:ref:`bg_homfractal`              Homogeneous distribution in fractal dimensions
==============================  =============================================================================

Phenomenological models
...........................................

This category comprises phenomenological models that represent various mathematical functions that are intended to mimick the background decay, without reference to a particular spatial distribution of spins. The parameters of these models do no have direct physical meaning.

.. rst-class:: func-list

==============================  =============================================================================
Model function                                  Description
==============================  =============================================================================
:ref:`bg_exp`                    Exponential decay
:ref:`bg_strexp`                 Stretched exponential decay
:ref:`bg_sumstrexp`              Sum of two stretched exponential decays
:ref:`bg_prodstrexp`             Product of two stretched exponential decays
:ref:`bg_poly1`                  Linear polynomial
:ref:`bg_poly2`                  Quadratic polynomial
:ref:`bg_poly3`                  Cubic polynomial
==============================  =============================================================================


.. _modelsref_ex:


Experiment models
------------------------------------

.. rst-class:: func-list

============================== =============================================================================
Model function                               Description
============================== =============================================================================
:ref:`ex_4pdeer`                  4-pulse DEER, 3-pulse DEER (single modulated pathway)
:ref:`ex_ovl4pdeer`               4-pulse DEER (including 2+1 pathway) 
:ref:`ex_5pdeer`                  5-pulse DEER
:ref:`ex_7pdeer`                  7-pulse DEER
:ref:`ex_ridme1`                  RIDME with one harmonic pathway (spin S=1/2)
:ref:`ex_ridme3`                  RIDME with three harmonic pathways (spin S=3/2)
:ref:`ex_ridme5`                  RIDME with five harmonic pathways (spin S=5/2)
:ref:`ex_ridme7`                  RIDME with seven harmonic pathways (spin S=7/2)
============================== =============================================================================



.. toctree::
    :maxdepth: 0
    :hidden:
    :includehidden:
    :glob:
    :caption: Models

    ./models/dd_*
    ./models/bg_*
    ./models/ex_*