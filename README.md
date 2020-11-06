<p align="center">
<img src="https://raw.githubusercontent.com/JeschkeLab/DeerLab/main/docsrc/source/_static/logo_dark.png" alt="DeerLab Logo" width="40%"></img>
</p>
</div>

### About
DeerLab is a Python package for the analysis of dipolar EPR (electron paramagnetic resonance) spectroscopy data. Dipolar EPR spectroscopy techniques include DEER (double electron-electron resonance), RIDME (relaxation-induced dipolar modulation enhancement), and others. The documentation can be found [here](https://jeschkelab.github.io/DeerLab/index.html).

DeerLab consists of a collection of functions for modelling, data processing, and least-squares fitting. They can be combined in scripts to build custom data analysis workflows. DeerLab supports both classes of distance distribution models: non-parametric (Tikhonov regularization and related) and parametric (multi-Gaussians etc). It also provides a selection of background and experiment models.

The early versions of DeerLab (up to version 0.9) are written in MATLAB. The old MATLAB codebase is archived and can be found [here](https://github.com/JeschkeLab/DeerLab-Matlab).

### Requirements

DeerLab is available for Windows, Mac and Linux systems and requires **Python 3.6**, **3.7**, or **3.8**.

All additional dependencies are automatically downloaded and installed during the setup.
 
### Setup

A pre-built distribution can be installed using `pip`.

First, ensure that `pip` is up-to-date. From a terminal (preferably with administrative privileges) use the following command:

    python -m pip install --upgrade pip

Next, install DeerLab with

    python -m pip install deerlab

More details on the installation of DeerLab can be found [here](https://jeschkelab.github.io/DeerLab/installation.html).

### Citation

When you use DeerLab in your work, please cite the following publication:

<div style="height:6%px; width:100%; display:flex; flex-wrap:wrap; align-items:center; border-style:solid">
    <div style="margin-left:2%; margin-bottom:2%; font-size:14px">
        <h3 style="font-size:110%">  DeerLab: a comprehensive software package for analyzing dipolar electron paramagnetic resonance spectroscopy data </h3> 
        Luis Fábregas Ibáñez, Gunnar Jeschke, Stefan Stoll <br>
        Magn. Reson., 1, 209–224, 2020 <br>
        <a href="https://doi.org/10.5194/mr-1-209-2020"> doi.org/10.5194/mr-1-209-2020</a>
    </div>
</div>

### License

The DeerLab toolbox is licensed under the [MIT License](LICENSE).

Copyright (c) 2019-2020: Luis Fábregas Ibáñez, Stefan Stoll, Gunnar Jeschke, and [other contributors](https://github.com/JeschkeLab/DeerLab/contributors).
