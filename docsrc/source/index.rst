DeerAnalysis2
=========================================

Overview
---------------------------------

This software package is a MATLAB toolbox for the analysis of data from DEER (double electron-electron resonance) spectroscopy and similar dipolar spectroscopy techniques (DQC, RIDME, SIFTER).

It consists of a collection of functions (the API) that provides a series of functions that perform single processing or fitting tasks. They can be combined in scripts to generate custom analysis workflows.

To model distance distributions, DeerAnalysis2 supports two types of model classes and associated workflows: parameter-free models (as used in Tikhonov regularization) as well as a series of parameterized models (mutli-Gaussians etc). It also provides a selection of background models. There are API functions for generating synthetic datasets as well as for fitting and analyzing experimental data sets.


Contents
==================

.. toctree::
    :maxdepth: 1

    ./tutorial
    ./functions
    ./examples
    ./theory


.. Indices and tables
.. ==================

.. * :ref:`genindex`
.. * :ref:`modindex`
.. * :ref:`search`
