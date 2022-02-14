.. _modelsref:


Models
======================

DeerLab provides a wide selection of built-in models. These are instances of the :ref:`Model` object class 
and can be manipulated and expanded as described in the :ref:`Modelling guide <modelling_guide>`.

.. currentmodule:: deerlab

.. _modelsref_dd:

These are models particularly oriented towards the modelling of dipolar EPR spectroscopy data. 

Distance distribution models 
****************************

These models provide parametrizations of the interspin distance distributions. 

.. rubric:: Basis functions


This group of distribution models is based on functions that can serve as basis functions for modeling more general distributions.

.. autosummary:: 
  :toctree: _autosummary
  :template: custom_data_template.rst
  :nosignatures:

  dd_gauss 
  dd_gauss2 
  dd_gauss3  
  dd_rice 
  dd_rice2 
  dd_rice3
  dd_gengauss
  dd_skewgauss
  dd_cos

.. rubric:: Disordered physical models

This group of distribution models represents spin labels attached to disordered macromolecular structures.

.. autosummary:: 
  :toctree: _autosummary
  :template: custom_data_template.rst
  :nosignatures:

  dd_wormchain
  dd_wormgauss
  dd_randcoil



.. rubric:: Sphere/shell models

This group of distribution models represents spin labels in simple partitions of 3D space.

.. autosummary:: 
  :toctree: _autosummary
  :template: custom_data_template.rst
  :nosignatures:

  dd_sphere
  dd_spheresurf
  dd_spherepoint
  dd_shell 
  dd_shellshell
  dd_shellsphere
  dd_shellvoidshell
  dd_shellvoidsphere



.. rubric:: Toy models

This group contains distribution models that have absolutely no physical relevance. They are useful for testing numerical algorithms.

.. autosummary:: 
  :toctree: _autosummary
  :template: custom_data_template.rst
  :nosignatures:

  dd_circle
  dd_triangle
  dd_uniform



.. _modelsref_bg:

Background models
******************

These models provide parametrizations of the signal contribution arising from intermolecular dipolar interactions.

.. rubric:: Physical models

This category involves models that describe particular distributions of spin labels in space. These models depend on physical parameters such as spin concentration, exclusion distances, and dimensionality.

.. autosummary:: 
  :toctree: _autosummary
  :template: custom_data_template.rst
  :nosignatures:

  bg_hom3d
  bg_hom3dex
  bg_homfractal


Additionally, the following models describe the time-dependent phase-shifts arising from these particular distributions of spin labels in space. 

.. autosummary:: 
  :toctree: _autosummary
  :template: custom_data_template.rst
  :nosignatures:

  bg_hom3d_phase
  bg_hom3dex_phase
  bg_homfractal_phase

.. rubric:: Phenomenological models

This category comprises phenomenological models that represent various mathematical functions that are intended to mimick the background decay, without reference to a particular spatial distribution of spins. The parameters of these models do no have direct physical meaning.

.. autosummary:: 
  :toctree: _autosummary
  :template: custom_data_template.rst
  :nosignatures:

  bg_exp
  bg_strexp
  bg_sumstrexp
  bg_prodstrexp
  bg_poly1
  bg_poly2
  bg_poly3


.. _modelsref_ex:


Experiment models
*******************

These generator functions model different dipolar EPR experiments, and provide the means
to refine the dipolar signal models by introducing experimental information. 

.. autosummary:: 
  :toctree: _autosummary
  :template: custom_function_template.rst
  :nosignatures:

  ex_3pdeer
  ex_4pdeer
  ex_rev5pdeer
  ex_fwd5pdeer
  ex_sifter
  ex_ridme