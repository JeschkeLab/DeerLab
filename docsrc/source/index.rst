DeerLab-development Documentation
=====================================================

.. warning:: DeerLab is currently in the pre-release stage (version numbers 0.x) and under active development. Major changes are likely before the first stable version (1.0) is released.

DeerLab is a MATLAB toolbox for the analysis of data from dipolar EPR (electron paramagnetic resonance) spectroscopy. Dipolar EPR spectroscopy techniques include DEER (double electron-electron resonance), RIDME (relaxation-induced dipolar modulation enhancement), and others.

DeerLab consists of a collection of functions that perform individual tasks. There are functions for generating synthetic signals as well as for processing and fitting experimental data.

DeerLab provides two types of distance distribution models: non-parametric models (as used in Tikhonov regularization) as well as a series of parametric distance distribution models (multi-Gaussians, multi-Rice, worm-like chain, etc.). It also provides a selection of background models (homogeneous spin distribution in 3 or fractal dimensions, excluded-volume effects). 

A publication about DeerLab is in preparation. When you use DeerLab in your work, please cite 

    Fábregas Ibáñez, L., Jeschke, G., and Stoll, S.: *DeerLab: A comprehensive toolbox for analyzing dipolar EPR spectroscopy data*, Magn. Reson. Discuss., https://doi.org/10.5194/mr-2020-13, in review, 2020

Please check back frequently for updated publication information.

.. toctree::
    :hidden:
    :caption: Getting DeerLab
    :maxdepth: 0

    ./download
    ./installation

.. toctree::
    :hidden:
    :caption: User Guide
    :maxdepth: 0

    ./firstscript
    ./basics
    ./preprocess
    ./fitting
    ./uncertainty
    ./examples

.. toctree::
    :hidden:
    :caption: Reference
    :maxdepth: 0

    ./functions
    ./modelsref
    ./theory

.. toctree::
    :hidden:
    :caption: Others
    :maxdepth: 1

    ./funding


.. Indices and tables
.. ==================

.. * :ref:`genindex`
.. * :ref:`modindex`
.. * :ref:`search`
