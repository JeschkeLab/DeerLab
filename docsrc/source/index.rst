DeerLab-development Documentation
=====================================================
       
Welcome to the documentation for DeerLab.

.. warning:: DeerLab is currently in the pre-release stage and under active development. Major changes are likely before the first stable version is released.

DeerLab is a MATLAB toolbox for the analysis of data from dipolar EPR (electron paramagnetic resonance) spectroscopy. Dipolar EPR spectroscopy techniques include DEER (double electron-electron resonance), RIDME (relaxation-induced dipolar modulation enhancement), and others.

DeerLab consists of a collection of functions that perform single processing or data fitting tasks. They can be combined in scripts to generate custom analysis workflows.

DeerLab provides two types of distance distribution models: parameter-free models (as used in Tikhonov regularization) as well as a series of parametrized distance distribution models (multi-Gaussians, worm-like chain, etc.). It also provides a selection of background models . There are functions for generating synthetic datasets as well as for fitting and analyzing experimental data sets.

.. toctree::
    :hidden:
    :caption: Getting DeerLab
    :maxdepth: 1

    ./download
    ./installation

.. toctree::
    :hidden:
    :caption: User Guide
    :maxdepth: 1

    ./firstscript
    ./basics
    ./models
    ./preprocess
    ./fitting
    ./uncertainty

.. toctree::
    :hidden:
    :caption: Examples
    :maxdepth: 1

    ./examples1
    ./examples2


.. toctree::
    :hidden:
    :caption: Reference
    :maxdepth: 1

    ./functions
    ./modelsref
    ./theory


.. Indices and tables
.. ==================

.. * :ref:`genindex`
.. * :ref:`modindex`
.. * :ref:`search`
