DeerLab-development Documentation
=====================================================
       
Welcome to the documentation for DeerLab.

.. note:: DeerLab is currently in its pre-release stage. Major changes (change of function names and calling syntax) are possible until the first stable version is released.

DeerLab is a MATLAB toolbox for the analysis of data from dipolar EPR (electron paramagnetic resonance) spectroscopy. This includes techniques such as DEER (double electron-electron resonance), RIDME (relaxation-induced dipolar modulation enhancement), and others.

DeerLab consists of a collection of functions that perform single processing or fitting tasks. They can be combined in scripts to generate custom analysis workflows. Many different models for time-domain signals and for distance distributions are available.

To model distance distributions, DeerLab supports two types of model classes and associated workflows: parameter-free models (as used in Tikhonov regularization) as well as a series of :doc:`parametrized distance distribution models <ddmodels>` (multi-Gaussians, worm-like chain, etc.). It also provides a selection of :doc:`background models <bgmodels>`. There are functions for generating synthetic datasets as well as for fitting and analyzing experimental data sets.

.. toctree::
    :hidden:
    :maxdepth: 1


    Home<self>

.. toctree::
    :hidden:
    :caption: User guide
    :maxdepth: 1

    ./download
    ./installation
    ./examples

.. toctree::
    :hidden:
    :caption: Reference
    :maxdepth: 1

    ./functions
    ./ddmodels
    ./bgmodels
    ./exmodels


.. Indices and tables
.. ==================

.. * :ref:`genindex`
.. * :ref:`modindex`
.. * :ref:`search`
