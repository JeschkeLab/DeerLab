Welcome to...
=========================================


.. image:: logo_dark.png
    :width: 650px
    :align: center

.. raw:: html

	<p align="center" style="top-padding: 100px">
	 <img src="https://img.shields.io/github/issues-raw/luisfabib/deerlab?style=flat"></img>
	 <img src="https://img.shields.io/badge/MATLAB-R2016b--R2019b-brightgreen?style=flat"></img>
	 <img src="https://img.shields.io/github/downloads/luisfabib/deerlab/total?style=flat"></img>
	 <img src="https://github.com/luisfabib/deerlab/workflows/Webpage%20update/badge.svg?style=flat-square"></img>
	 <img src="https://img.shields.io/endpoint?url=https%3A%2F%2Fdeershields.s3.eu-central-1.amazonaws.com%2Fcoverage_badge.json"></img>
	 <img src="https://img.shields.io/endpoint?url=https%3A%2F%2Fdeershields.s3.eu-central-1.amazonaws.com%2Ftestsuite_badge.json"></img>
	</p>

       
-------

.. raw:: html
     :file: ./_static/downloadbutton.html
       
-------

.. important:: DeerLab is currently in its pre-release stage. Major changes (change of function names and calling syntax) are possible until the first stable version is released.

DeerLab is a MATLAB toolbox for the analysis of data from dipolar EPR (electron paramagnetic resonance) spectroscopy. This includes techniques such as DEER (double electron-electron resonance), RIDME (relaxation-induced dipolar modulation enhancement), and others.

DeerLab consists of a collection of functions that perform single processing or fitting tasks. They can be combined in scripts to generate custom analysis workflows. Many different models for time-domain signals and for distance distributions are available.

To model distance distributions, DeerLab supports two types of model classes and associated workflows: parameter-free models (as used in Tikhonov regularization) as well as a series of parameterized models (mutli-Gaussians, worm-like chain, etc.). It also provides a selection of background models. There are functions for generating synthetic datasets as well as for fitting and analyzing experimental data sets.


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
    ./models


.. Indices and tables
.. ==================

.. * :ref:`genindex`
.. * :ref:`modindex`
.. * :ref:`search`
