# DeerAnalysis2

<a href="https://github.com/luisfabib/DeerAnalysis2/issues"><img alt="GitHub issues" src="https://img.shields.io/github/issues/luisfabib/DeerAnalysis2"></a>
<img src="https://img.shields.io/github/license/mashape/apistatus.svg"></img>
<img alt="GitHub top language" src="https://img.shields.io/github/languages/top/luisfabib/DeerAnalysis2">
<img alt="GitHub All Releases" src="https://img.shields.io/github/downloads/luisfabib/DeerAnalysis/total">

This software package is a MATLAB toolbox for the analysis of data from DEER (double electron-electron resonance) spectroscopy and similar dipolar spectroscopy techniques (DQC, RIDME, SIFTER).

It consists of a collection of functions (the API) that provides a series of functions that perform single processing or fitting tasks. They can be combined in scripts to generate custom analysis workflows.

To model distance distributions, DeerAnalysis2 supports two types of model classes and associated workflows: parameter-free models (as used in Tikhonov regularization) as well as a series of parameterized models (mutli-Gaussians etc). It also provides a selection of background models. There are API functions for generating synthetic datasets as well as for fitting and analyzing experimental data sets.
