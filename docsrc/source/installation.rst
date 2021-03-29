.. _installation:

Installation
=========================

Requirements
------------

To install DeerLab, first install Python on your computer. Python can be downloaded from the `official Python distribution <https://www.python.org/>`_. There are
many online tutorials to guide you through the installation and setup (see `here <https://realpython.com/installing-python/>`_ for example). Make sure you install
one of the Python versions compatible with DeerLab, either **Python 3.6**, **3.7**, **3.8**, or **3.9**  

For Windows systems it is important to ensure that the **Install launcher for all users (recommended)** and  the **Add Python 3.x to PATH** checkboxes at the bottom are checked. To test if python has been successfully  installed, open a terminal window and run the command::

	python

which should display the installed Python version and launch the Python command line. To exit, use the ``exit()`` command.

Installing pre-built DeerLab
-----------------------------

A pre-built distribution can be installed using `pip <https://pip.pypa.io/en/stable/installing/>`_.

First, it is important to ensure that ``pip`` is up-to-date. From a terminal (preferably with administrative privileges) use the following command::

		python -m pip install --upgrade pip

Next, install DeerLab and all its dependencies via::

		python -m pip install deerlab

DeerLab installs the following packages:

	* `matplotlib <https://matplotlib.org/>`_ - A comprehensive library for creating data visualizations with Python
	* `memoization <https://pypi.org/project/memoization/>`_ - A powerful caching library for Python
	* `pytest <https://docs.pytest.org/en/stable/>`_ - A Python testing framework
	* `cvxopt <https://cvxopt.org/index.html>`_ - Free software package for convex optimization
	* `numpy <https://numpy.org/>`_ -  Base N-dimensional array package 
	* `scipy <https://www.scipy.org/>`_ - Fundamental library for scientific computing
	* `joblib <https://joblib.readthedocs.io/en/latest/>`_ - Library lightweight pipelining and parallelization.

The installed numerical packages (numpy, scipy, cvxopt) are linked against different BLAS libraries depending on the OS:

	* Windows: linked against the Intel Matrix Kernel Library (MKL)
	* Linux: linked against OpenBLAS
	* Mac: linked against BLAS/LAPACK from the Accelerate framework

If an error occurs during or after the installation, please consult `this section <./installation_failed.html>`_ for a possible solution.

Upgrading to the latest version 
--------------------------------
To upgrade an existing DeerLab installation to the latest released version, use the following command from a terminal:: 

		python -m pip install --upgrade deerlab


Installing specific versions
-----------------------------

Any DeerLab version released after v0.10.0 can be installed via pip using the following command matching the x.y.z to the desired version::

		python -m pip install deerlab==x.y.z


DeerLab version prior to 0.10 are written in MATLAB and are still available from an `archived repository <https://github.com/JeschkeLab/DeerLab-Matlab>`_. 
Download and installation instruction for the MATLAB environment are provided in the released documentation. MATLAB releases have been deprecated and no further support is provided.


Installing from source
----------------------

To install DeerLab from the source, downloaded or clone the source code from the `DeerLab repository <https://github.com/JeschkeLab/DeerLab>`_.
DeerLab and its dependencies can be installed by running the following command on a terminal window to install DeerLab as a static Python package::

		python -m setup.py install


For developers, in order to install DeerLab but be able to frequently edit the code without having to re-install the package every time, use the command::

		python -m setup.py develop


Any changes made to the source code will then immediate effect.
