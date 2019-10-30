

.. image:: logo_dark.png
    :width: 800px
    :align: center


.. image:: https://img.shields.io/github/issues-raw/luisfabib/DeerAnalysis2?style=flat-square
.. image:: https://img.shields.io/badge/MATLAB-%3C2019a-brightgreen?style=flat-square
.. image:: https://img.shields.io/github/downloads/luisfabib/DeerAnalysis2/total?style=flat-square
.. image:: https://img.shields.io/badge/license-MIT-blue?style=flat-square


       
-------

.. raw:: html
     :file: ./_static/downloadbutton.html
	
       
-------



This software package is a MATLAB toolbox for the analysis of data from DEER (double electron-electron resonance) spectroscopy and similar dipolar spectroscopy techniques (DQC, RIDME, SIFTER).




It consists of a collection of functions (the API) that provides a series of functions that perform single processing or fitting tasks. They can be combined in scripts to generate custom analysis workflows.

To model distance distributions, DeerAnalysis2 supports two types of model classes and associated workflows: parameter-free models (as used in Tikhonov regularization) as well as a series of parameterized models (mutli-Gaussians etc). It also provides a selection of background models. There are API functions for generating synthetic datasets as well as for fitting and analyzing experimental data sets.


.. toctree::
    :caption: First Steps
    :maxdepth: 1

    ./download
    ./tutorial

.. toctree::
    :caption: API Documentation
    :maxdepth: 1

    ./functions
    ./models
    ./examples


.. toctree::
    :caption: Theory
    :maxdepth: 1

    ./theory


.. Indices and tables
.. ==================

.. * :ref:`genindex`
.. * :ref:`modindex`
.. * :ref:`search`
