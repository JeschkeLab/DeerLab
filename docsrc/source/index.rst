DeerLab |version|
=====================================================

.. warning:: DeerLab is currently in the pre-release stage (version numbers 0.x.x) and under active development. Major changes are likely before the first stable version (1.0.0) is released.

DeerLab is a free software package for the analysis of dipolar EPR (electron paramagnetic resonance) spectroscopy data based on the Python programming language. Dipolar EPR spectroscopy techniques include DEER (double electron-electron resonance), RIDME (relaxation-induced dipolar modulation enhancement), and others.

DeerLab consists of a collection of functions that perform modelling, processing or fitting tasks. They can be combined in scripts to build custom data analysis workflows. To model distance distributions, DeerLab supports two types of model classes and associated workflows: parameter-free models (as used in Tikhonov regularization) as well as a series of parameterized models (multi-Gaussians etc). It also provides a selection of background and experiment models. 

-------

Documentation
    To get started read the :ref:`beginners_guide` tutorial. 

    Looking for a specific topic? Look at the :ref:`selection of examples <examples>` or at :ref:`reference <functions>` for specific DeerLab functions. 

-------

What's new?
    Check the :ref:`changelog` for the changes introduced in the latest releases.

-------

Citing DeerLab
    When you use DeerLab in your work, please cite the following publication:

    .. raw:: html 

        <div style="height:100px; width:100%; display:flex; flex-wrap:wrap; align-items:center;">
            <div style="height:100%;">
                <img src='https://www.magnetic-resonance-ampere.net/graphic_journal_cover_normal.png', style="margin-top:1%; margin-bottom:1%;height:96%;">
            </div>
            <div style="margin-left:1%; font-size:14px">
                <h3 style="font-size:110%">  DeerLab: a comprehensive software package for analyzing dipolar electron paramagnetic resonance spectroscopy data </h3> 
                Luis Fábregas Ibáñez, Gunnar Jeschke, Stefan Stoll <br>
                Magn. Reson., 1, 209–224, 2020 <br>
                <a href="https://doi.org/10.5194/mr-1-209-2020"> doi.org/10.5194/mr-1-209-2020</a>
            </div>
        </div>


-------

Contributing
    .. raw:: html 

        <div style="height:50px; width:100%; display:flex; flex-wrap:wrap; align-items:center;">
            <div style="height:100%;">
                <img src='https://upload.wikimedia.org/wikipedia/commons/thumb/9/91/Octicons-mark-github.svg/1200px-Octicons-mark-github.svg.png', style="margin-top:1%; margin-bottom:1%;height:96%;">
            </div>
            <div style="margin-left:1%; font-size:14px">
                DeerLab is hosted on <a href="https://github.com/JeschkeLab/DeerLab"> Github </a>:
                <ul style="font-size:14px">
                <li> Bugs and feature requests can be reported on the <a href="https://github.com/JeschkeLab/DeerLab/issues"> issue tracker </a>.
                <li> Any contributions in the form of <a href="https://github.com/JeschkeLab/DeerLab/pulls"> pull requests</a>  are very welcome.
                </ul>
            </div>
        </div>


.. include:: changelog.rst
    :start-after: Version

.. toctree::
    :hidden:
    :caption: User Guide
    :maxdepth: 0

    ./installation
    ./basics
    ./beginners_guide
    ./uncertainty
    ./examples

.. toctree::
    :hidden:
    :caption: Reference
    :maxdepth: 0

    ./functions
    ./modelsref
    ./theory

.. toctree::
    :hidden:
    :caption: Others
    :maxdepth: 1

    ./license
    ./funding
    ./changelog

.. Indices and tables
.. ==================

.. * :ref:`genindex`
.. * :ref:`modindex`
.. * :ref:`search`
