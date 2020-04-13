Distribution models
=========================================

In DeerLab, distance distributions can be represented either as parameter-free models or as parametric models.

Parameter-free distance distributions
------------------------------------------------

.. rst-class:: coderef

A parameter-free distance distribution is a vector of population densities ``P`` defined over a vector of distances ``r``. Outside the range of ``r``, the distance distribution is considered zero, i.e. it is truncated to the range of ``r``. Such parameter-free distance distributions are returned by the fit functions :ref:`fitregmodel` and :ref:`obir`. A ``P`` can be converted to a time-domain signal using :ref:`dipolarsignal` or :ref:`dipolarkernel`.

Parameter-free distance distributions are preferred over parametric distance distributions, since they make fewer assumptions about the distribution and are more flexible. They introduce less bias.

Parametric distance distributions
------------------------------------------------

.. rst-class:: coderef

DeerLab provide a series of parametric distance distributions models. These can be used with several functions, including :ref:`fitparamodel` and :ref:`selectmodel`. All parametric distance distribution models are functions start with the prefix ``dd_``. They take a distance vector ``r`` and a parameter vector ``param`` as inputs and return the distance distribution as a vector ``P``. Here is an example:

.. code-block:: matlab

        P = dd_rice(r,param)


DeerLab's parametric distance distribution models can be classified into several groups. There are models that are based on basis functions and linear combinations thereof. Gaussian and 3D-Rice functions are the most important ones.

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
============================= ===================================================

A second group of distance distribution models represent three-dimensional disordered segmented objects such as proteins and other polymers.

.. rst-class:: func-list

=========================== ===================================================
  Model function                                       Description
=========================== ===================================================
:ref:`dd_wormchain`           Worm-like chain
:ref:`dd_wormgauss`           Worm-like chain with Gaussian convolution
:ref:`dd_randcoil`            Random coil
=========================== ===================================================


A third group provides models for distributions of labels in simple confined spaces such as spheres and spherical shells.

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


Finally, DeerLab provides a series of distance distribution models that represent simple geometric shapes. These models are worthless for practical purposes, but are useful as toy distributions for methods development.

.. rst-class:: func-list

============================== =============================================================================
  Model function                                       Description
============================== =============================================================================
:ref:`dd_circle`                Semi-circle distribution
:ref:`dd_triangle`              Triangle distribution
:ref:`dd_uniform`               Uniform distribution
============================== =============================================================================


.. toctree::
    :maxdepth: 0
    :hidden:
    :glob:
    :caption: Background models

    ./models/dd_*
