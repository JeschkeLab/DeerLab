# DeerLab

<p align="center">
<img src="./docsrc/source/_static/logo_dark.png" alt="DeerLab Logo" width="60%"></img>
</p>
</div>

### About
DeerLab is a free software package for the analysis of dipolar EPR (electron paramagnetic resonance) spectroscopy data based on the Python programming language. Dipolar EPR spectroscopy techniques include DEER (double electron-electron resonance), RIDME (relaxation-induced dipolar modulation enhancement), and others. The main homepage can be found at [jeschkelab.github.io/DeerLab](https://jeschkelab.github.io/DeerLab/). This is the GitHub repository of the DeerLab source code, including instructions for compiling and installing DeerLab.

DeerLab consists of a collection of functions that perform modelling, processing or fitting tasks. They can be combined in scripts to build custom data analysis workflows. To model distance distributions, DeerLab supports two types of model classes and associated workflows: parameter-free models (as used in Tikhonov regularization) as well as a series of parameterized models (mutli-Gaussians etc). It also provides a selection of background and experiment models.

### Requirements

DeerLab requires **Python 3.6-3.8** and is available for Windows, Mac and Linux systems.

All additional dependencies will be automatically downloaded and installed during the setup.
 
### Setup

A pre-built distribution can be installed using pip, which is already installed in Python 3.4 or newer versions.

First it is important to ensure that pip is up-to-date. From a terminal (preferably with administrative privileges) use the following command::

    python -m pip install --upgrade pip

Next, DeerLab and all its dependencies can be installed via::

    python -m pip install deerlab

More details on the installation of DeerLab can be found [here](https://jeschkelab.github.io/DeerLab/develop/installation.html).

### Citation

A publication about DeerLab is [available here](https://doi.org/10.5194/mr-2020-13). When you use DeerLab in your work, please cite 

> Fábregas Ibáñez, L., Jeschke, G., and Stoll, S.: DeerLab: A comprehensive toolbox for analyzing dipolar EPR spectroscopy data, Magn. Reson. Discuss., https://doi.org/10.5194/mr-2020-13, 2020

Please check back frequently for updated publication information.

### License

The DeerLab toolbox is licensed under the MIT License.

Copyright (c) 2019-2020: Luis Fábregas Ibáñez, Stefan Stoll, Gunnar Jeschke, and [other contributors](https://github.com/JeschkeLab/DeerLab/contributors).
