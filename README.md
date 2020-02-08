
<p align="center">
<img src="./docsrc/source/logo_dark.png" alt="DeerLab Logo" width="60%"></img>
</p>
</div>

<p align="center">
  <img src="https://img.shields.io/github/issues-raw/luisfabib/deerlab?style=flat"></img>
  <img src="https://img.shields.io/badge/MATLAB-R2016b--R2019b-brightgreen?style=flat"></img>
  <img src="https://img.shields.io/github/downloads/luisfabib/deerlab/total?style=flat"></img>
  <img src="https://github.com/luisfabib/deerlab/workflows/Webpage%20update/badge.svg?style=flat-square"></img>
  <img src="https://img.shields.io/endpoint?url=https%3A%2F%2Fdeershields.s3.eu-central-1.amazonaws.com%2Fcoverage_badge.json"></img>
  <img src="https://img.shields.io/endpoint?url=https%3A%2F%2Fdeershields.s3.eu-central-1.amazonaws.com%2Ftestsuite_badge.json"></img>
</p>

### About
The DeerLab software package is a MATLAB toolbox for the analysis of data from DEER (double electron-electron resonance) spectroscopy and similar dipolar EPR spectroscopy techniques (RIDME, DQC, SIFTER). The main homepage can be found at www.deeranalysis.org. This is the GitHub repository of the DeerLab source code, including instructions for compiling and installing DeerLab.

It consists of a collection of functions that perform processing or fitting tasks. They can be combined in scripts to build custom data analysis workflows.

To model distance distributions, DeerLab supports two types of model classes and associated workflows: parameter-free models (as used in Tikhonov regularization) as well as a series of parameterized models (mutli-Gaussians etc). It also provides a selection of background models. There are functions for generating synthetic datasets as well as for fitting and analyzing experimental data sets.

### Requirements
DeerLab requires the following products:

  * MATLAB (R2016b or newer) (see <https://ch.mathworks.com/products/matlab.html>)
 
 Optional functionality may require the following products:
 
  * Optimization Toolbox (see <https://ch.mathworks.com/products/optimization.html>)

### Setup

In order for MATLAB to access the DeerLab API functions, the path to the DeerLab installation folder must be set in MATLAB.

**Option 1:** Add DeerLab path via MATLAB's IDE

1) On the ``Home`` tab, in the ``Environment`` section, click ``Set Path``. 

2) Click ``Add with Subfolders...`` and select the ``DeerLab\functions`` directory. 

3) Click ``Save`` to save the current MATLAB search path and exit via ``Close``.

**Option2:**  Add DeerLab path at startup

1) Open (or create) the ``startup.m`` file in the default ``\MATLAB`` directory.

2) Add the following lines of code:

       addpath('mypath/DeerLab/functions')
       addpath('mypath/DeerLab/functions/models')

3) Save ``startup.m`` and restart MATLAB.

### License

The DeerLab toolbox is licensed under the MIT License. The complete toolbox consists of the functions ([functions/](https://github.com/luisfabib/deerlab/tree/master/functions)), documentation source ([docsrc/](https://github.com/luisfabib/deerlab/tree/master/docsrc)), tutorial scripts ([tutorials/](https://github.com/luisfabib/deerlab/tree/master/tutorials)), test suite ([tests/](https://github.com/luisfabib/deerlab/tree/master/tests)), and pipeline scripts ([.github/workflows](https://github.com/luisfabib/deerlab/tree/master/.github/workflows)). See below for exceptions.

DeerLab includes code from the following projects, which have their own licenses:
- [datahash.m](https://www.mathworks.com/matlabcentral/fileexchange/31272-datahash) (Hash-key generator by Jan Simon) [BSD] 
- [fresnelS.m, fresnelC.m](https://www.mathworks.com/matlabcentral/fileexchange/28765-fresnels-and-fresnelc) (Efficient and accurate Fresnel integrals by John D'Errico) [BSD]
- [fminsearchcon.m](https://www.mathworks.com/matlabcentral/fileexchange/8277-fminsearchbnd-fminsearchcon) (Bound constrained optimization using fminsearch by John D'Errico) [BSD]
- [nlsqbnd.m](https://www.mathworks.com/matlabcentral/fileexchange/8277-fminsearchbnd-fminsearchcon) (Non-linear least squares solver with box constraints by Alain Barraud) [BSD]

Copyright (c) 2019-2020: Luis Fabregas, Stefan Stoll, Gunnar Jeschke, and [other contributors](https://github.com/luisfabib/deerlab/contributors).
