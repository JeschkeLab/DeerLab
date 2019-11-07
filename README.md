
<a name="logo"/>
<img src="./docsrc/source/logo_dark.png" alt="DeerAnalysis Logo" width="700" height="150"></img>
</a>
</div>

![GitHub issues](https://img.shields.io/github/issues-raw/luisfabib/DeerAnalysis2?style=flat)
<img src="https://img.shields.io/badge/MATLAB-R2016b--R2019b-brightgreen?style=flat"></img>
![GitHub All Releases](https://img.shields.io/github/downloads/luisfabib/DeerAnalysis2/total?style=flat)
![GItHub Actions](https://github.com/luisfabib/DeerAnalysis2/workflows/Webpage%20update/badge.svg?style=flat-square)

This DeerAnalysis2 software package is a MATLAB toolbox for the analysis of data from DEER (double electron-electron resonance) spectroscopy and similar dipolar spectroscopy techniques (DQC, RIDME, SIFTER). The main homepage for Julia can be found at www.deeranalysis.org. This is the GitHub repository of the DeerAnalysis2 source code, including instructions for compiling and installing DeerAnalysis.

It consists of a collection of functions (the API) that provides a series of functions that perform single processing or fitting tasks. They can be combined in scripts to generate custom analysis workflows.

To model distance distributions, DeerAnalysis2 supports two types of model classes and associated workflows: parameter-free models (as used in Tikhonov regularization) as well as a series of parameterized models (mutli-Gaussians etc). It also provides a selection of background models. There are API functions for generating synthetic datasets as well as for fitting and analyzing experimental data sets.

Requirements
---------------
The application programming interface (API) of DeerAnalysis requires the following products:

  * MATLAB (oldest version supported R2016b) (see <https://ch.mathworks.com/products/matlab.html>)
    
  * Optimization Toolbox (see <https://ch.mathworks.com/products/optimization.html>)

Setup
---------------
In order for MATLAB to access the DeerAnalysis API functions, the path to DeerAnalysis installation folder must be set.

**Option 1:** Add DeerAnalysis path via MATLAB's IDE

1) On the ``Home`` tab, in the ``Environment`` section, click ``Set Path``. 

2) Click ``Add with Subfolders...`` and select the ``DeerAnalysis\functions`` directory. 

3) Click ``Save`` to save the current MATLAB search path and exit via ``Close``.

**Option2:**  Add DeerAnalysis path at startup

1) Open (or create) the ``startup.m`` file in the default ``\MATLAB`` directory.

2) Add the following lines of code:

       addpath('mypath/DeerAnalysis/functions')
       addpath('mypath/DeerAnalysis/functions/models')

3) Save ``startup.m`` and restart MATLAB.


