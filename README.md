
<p align="center">
<img src="./docsrc/source/_static/logo_dark.png" alt="DeerLab Logo" width="60%"></img>
</p>
</div>

<p align="center">
  <img src="https://img.shields.io/github/issues-raw/JeschkeLab/DeerLab?style=flat"></img>
  <img src="https://img.shields.io/badge/MATLAB-R2017a--R2020a-brightgreen?style=flat"></img>
  <img src="https://img.shields.io/github/downloads/JeschkeLab/DeerLab/total?style=flat"></img>
  <img src="https://img.shields.io/endpoint?url=https%3A%2F%2Fdeershields.s3.eu-central-1.amazonaws.com%2Fcoverage_badge.json"></img>
  <img src="https://img.shields.io/endpoint?url=https%3A%2F%2Fdeershields.s3.eu-central-1.amazonaws.com%2Ftestsuite_badge.json"></img>
</p>

### About
The DeerLab software package is a MATLAB toolbox for the analysis of data from DEER (double electron-electron resonance) spectroscopy and similar dipolar EPR spectroscopy techniques (RIDME, DQC, SIFTER,...). The main homepage can be found at [jeschkelab.github.io/DeerLab](https://jeschkelab.github.io/DeerLab/). This is the GitHub repository of the DeerLab source code, including instructions for compiling and installing DeerLab.

DeerLab consists of a collection of functions that perform modelling, processing or fitting tasks. They can be combined in scripts to build custom data analysis workflows.

To model distance distributions, DeerLab supports two types of model classes and associated workflows: parameter-free models (as used in Tikhonov regularization) as well as a series of parameterized models (mutli-Gaussians etc). It also provides a selection of background and experiment models. There are functions for generating synthetic datasets as well as for fitting and analyzing experimental data sets.

If you want to get notified upon publication of a new DeerLab release, you just need to press the `Watch` button above and select `Releases only`. 

### Requirements
DeerLab requires the following products:

  * MATLAB (R2017a or newer) (see <https://www.mathworks.com/products/matlab.html>)
 
DeerLab will use the following product if installed:
 
  * Optimization Toolbox (see <https://www.mathworks.com/products/optimization.html>)

### Setup

In order for MATLAB to access the DeerLab functions, the path to the DeerLab installation folder must be set in MATLAB.

**Option 1:** Add DeerLab path via MATLAB's IDE

1) On the ``Home`` tab, in the ``Environment`` section, click ``Set Path``. 

2) Click ``Add with Subfolders...`` and select the ``DeerLab\functions`` directory. 

3) Click ``Save`` to save the current MATLAB search path and exit via ``Close``.

**Option2:**  Add DeerLab path at startup

1) Open (or create) the ``startup.m`` file in the default ``\MATLAB`` directory.

2) Add the following lines of code:

       addpath('mypath/DeerLab/functions')

3) Save ``startup.m`` and restart MATLAB.

### Citation

A publication about DeerLab is in preparation. When you use DeerLab in your work, please cite 

> L. Fábregas Ibáñez, G. Jeschke, S. Stoll, DeerLab: A comprehensive toolbox for analyzing dipolar EPR spectroscopy data, *in preparation*

Please check back frequently for updated publication information.

### License

The DeerLab toolbox is licensed under the MIT License. The complete toolbox consists of the functions ([functions/](https://github.com/JeschkeLab/DeerLab/tree/master/functions)), documentation source ([docsrc/](https://github.com/JeschkeLab/DeerLab/tree/master/docsrc)), tutorial scripts ([tutorials/](https://github.com/JeschkeLab/DeerLab/tree/master/tutorials)), test suite ([tests/](https://github.com/JeschkeLab/DeerLab/tree/master/tests)), and pipeline scripts ([.github/workflows](https://github.com/JeschkeLab/DeerLab/tree/master/.github/workflows)). See below for exceptions.

DeerLab includes code from the following projects, which have their own licenses:
- [datahash.m](https://www.mathworks.com/matlabcentral/fileexchange/31272-datahash) (Hash-key generator by Jan Simon) [BSD] 
- [fresnelS.m, fresnelC.m](https://www.mathworks.com/matlabcentral/fileexchange/28765-fresnels-and-fresnelc) (Efficient and accurate Fresnel integrals by John D'Errico) [BSD]
- [fminsearchcon.m](https://www.mathworks.com/matlabcentral/fileexchange/8277-fminsearchbnd-fminsearchcon) (Bound constrained optimization using fminsearch by John D'Errico) [BSD]
- [nlsqbnd.m](https://ch.mathworks.com/matlabcentral/fileexchange/23621-nlsqbnd) (Non-linear least squares solver with box constraints by Alain Barraud) [BSD]
- [golden.m](https://www.mathworks.com/matlabcentral/fileexchange/25919-golden-section-method-algorithm) (Golden Section method algorithm by Katarzyna Zarnowiec) [BSD]
- [jacobianest.m](https://www.mathworks.com/matlabcentral/fileexchange/13490-adaptive-robust-numerical-differentiation) (Adaptive Robust Numerical Differentiation by John D'Errico) [BSD]
- [kde.m](https://ch.mathworks.com/matlabcentral/fileexchange/14034-kernel-density-estimator) (Kernel Density Estimator by Zdravko Botev) [BSD]
- [LevenbergMarquardt.m, jacobiansimple.m](https://ch.mathworks.com/matlabcentral/fileexchange/53449-levenberg-marquardt-toolbox)(Levenberg-Marquardt & Jacobian toolbox by Alexander Dentler)[BSD]
- [fdcoeffF.m](https://faculty.washington.edu/rjl/fdmbook/matlab/fdcoeffF.m)(Fornberg's method for finite difference coefficients)

Copyright (c) 2019-2020: Luis Fábregas Ibáñez, Stefan Stoll, Gunnar Jeschke, and [other contributors](https://github.com/JeschkeLab/DeerLab/contributors).
