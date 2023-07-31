# DeerLab

[![https://jeschkelab.github.io/DeerLab/](https://img.shields.io/pypi/v/deerlab)](https://pypi.org/project/DeerLab/)
[![Website](https://img.shields.io/website?down_message=offline&label=Documentation&up_message=online&url=https%3A%2F%2Fjeschkelab.github.io%2FDeerLab%2Findex.html)](https://jeschkelab.github.io/DeerLab/)
[![PyPI - Python Version](https://img.shields.io/pypi/pyversions/deerlab)](https://www.python.org/downloads/)
![PyPI - Downloads](https://img.shields.io/pypi/dm/deerlab?color=brightgreen)

## About

DeerLab is a comprehensive free scientific software package for Python focused on modeling, penalized least-squares regression, and uncertainty quantification. 
It provides highly specialized on the analysis of dipolar EPR (electron paramagnetic resonance) spectroscopy data. Dipolar EPR spectroscopy techniques include DEER (double electron-electron resonance), RIDME (relaxation-induced dipolar modulation enhancement), and others. 

The documentation can be found [here](https://jeschkelab.github.io/DeerLab/index.html).

The early versions of DeerLab (up to version 0.9.2) are written in MATLAB. The old MATLAB codebase is archived and can be found [here](https://github.com/JeschkeLab/DeerLab-Matlab).

## Requirements

DeerLab is available for Windows, Mac and Linux systems and requires **Python 3.8**, **3.9**, **3.10**, or  **3.11**.

All additional dependencies are automatically downloaded and installed during the setup.
 
## Setup

A pre-built distribution can be installed from the PyPI repository using `pip`.

From a terminal (preferably with admin privileges) use the following command to install from PyPI:

    python -m pip install deerlab

More details on the installation and updating of DeerLab can be found [here](https://jeschkelab.github.io/DeerLab/installation.html).

## Citing DeerLab

When you use DeerLab in your work, please cite the following publication:

 **DeerLab: a comprehensive software package for analyzing dipolar electron paramagnetic resonance spectroscopy data** <br>
 Luis Fábregas Ibáñez, Gunnar Jeschke, Stefan Stoll <br>
 Magn. Reson., 1, 209–224, 2020 <br>
 <a href="https://doi.org/10.5194/mr-1-209-2020"> doi.org/10.5194/mr-1-209-2020</a>

Here is the citation in bibtex format:

```tex
@article{FabregasIbanez2020_DeerLab,
  title = {{DeerLab}: a comprehensive software package for analyzing dipolar electron paramagnetic resonance spectroscopy data},
  author = {Fábregas Ibáñez, Luis and Jeschke, Gunnar and Stoll, Stefan},
  journal = {Magnetic Resonance},
  year = {2020},
  volume = {1},
  number = {2},
  pages = {209--224},
  doi = {10.5194/mr-1-209-2020}
}
```

## License

DeerLab is licensed under the [MIT License](LICENSE).

Copyright © 2019-2023: Luis Fábregas Ibáñez, Stefan Stoll, Gunnar Jeschke, and [other contributors](https://github.com/JeschkeLab/DeerLab/contributors).
