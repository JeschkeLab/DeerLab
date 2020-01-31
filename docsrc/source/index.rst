Welcome to...
=========================================

.. image:: logo_dark.png
    :width: 800px
    :align: center

.. raw:: html

	<p align="center">
	 <img src="https://img.shields.io/github/issues-raw/luisfabib/DeerAnalysis2?style=flat"></img>
	 <img src="https://img.shields.io/badge/MATLAB-R2016b--R2019b-brightgreen?style=flat"></img>
	 <img src="https://img.shields.io/github/downloads/luisfabib/DeerAnalysis2/total?style=flat"></img>
	 <img src="https://github.com/luisfabib/DeerAnalysis2/workflows/Webpage%20update/badge.svg?style=flat-square"></img>
	 <img src="https://img.shields.io/endpoint?url=https%3A%2F%2Fdeershields.s3.eu-central-1.amazonaws.com%2Fcoverage_badge.json"></img>
	 <img src="https://img.shields.io/endpoint?url=https%3A%2F%2Fdeershields.s3.eu-central-1.amazonaws.com%2Ftestsuite_badge.json"></img>
	</p>


       
-------

.. raw:: html
     :file: ./_static/downloadbutton.html
       
-------

.. important:: DeerAnalysis is currently on its pre-release phase. Minor changes and some major changes can be expected in the short-term future until the first stable version is released. 

This software package is a MATLAB toolbox for the analysis of data from DEER (double electron-electron resonance) spectroscopy and similar dipolar spectroscopy techniques (DQC, RIDME, SIFTER).




It consists of a collection of functions (the API) that provides a series of functions that perform single processing or fitting tasks. They can be combined in scripts to generate custom analysis workflows.

To model distance distributions, DeerAnalysis2 supports two types of model classes and associated workflows: parameter-free models (as used in Tikhonov regularization) as well as a series of parameterized models (mutli-Gaussians etc). It also provides a selection of background models. There are API functions for generating synthetic datasets as well as for fitting and analyzing experimental data sets.


.. toctree::
    :hidden:
    :caption: First Steps
    :maxdepth: 1

    ./download
    ./installation

.. toctree::
    :hidden:
    :caption: API Documentation
    :maxdepth: 1

    ./functions
    ./models
    ./examples


.. toctree::
    :hidden:
    :caption: Theory
    :maxdepth: 1

    ./theory


.. Indices and tables
.. ==================

.. * :ref:`genindex`
.. * :ref:`modindex`
.. * :ref:`search`
