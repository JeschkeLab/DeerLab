
<a name="logo"/>
<img src="./docsrc/source/favicon.png" alt="DeerAnalysis Logo" width="100" height="100"></img>
</a>
</div>

# DeerAnalysis2
![GitHub issues](https://img.shields.io/github/issues-raw/luisfabib/DeerAnalysis2?style=flat-square)
<img src="https://img.shields.io/badge/MATLAB-%3C2019a-brightgreen?style=flat-square"></img>
![GitHub All Releases](https://img.shields.io/github/downloads/luisfabib/DeerAnalysis2/total?style=flat-square)
<img src="https://img.shields.io/badge/license-MIT-blue?style=flat-square"></img>
!![](https://github.com/luisfabib/DeerAnalysis2/workflows/Webpage%20update/badge.svg)

This software package is a MATLAB toolbox for the analysis of data from DEER (double electron-electron resonance) spectroscopy and similar dipolar spectroscopy techniques (DQC, RIDME, SIFTER).

It consists of a collection of functions (the API) that provides a series of functions that perform single processing or fitting tasks. They can be combined in scripts to generate custom analysis workflows.

To model distance distributions, DeerAnalysis2 supports two types of model classes and associated workflows: parameter-free models (as used in Tikhonov regularization) as well as a series of parameterized models (mutli-Gaussians etc). It also provides a selection of background models. There are API functions for generating synthetic datasets as well as for fitting and analyzing experimental data sets.
