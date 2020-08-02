Installation instructions
=========================

Requirements
------------
DeerLab requires one of the following versions of the Python interpreter
	
	* Python 3.6.x
	* Python 3.7.x
	* Python 3.8.x

which can be downloaded from the `official Python distribution <https://www.python.org/>`_.

For Windows systems it is imporant to ensure that the **Install launcher for all users (recommended)** and the **Add Python 3.7 to PATH** checkboxes at the bottom are checked. 

Installing pre-built DeerLab
-----------------------------
A pre-built distribution can be installed using `pip <https://pip.pypa.io/en/stable/installing/>`_, which is already installed in Python 3.4 or newer versions. 

First it is important to ensure that ``pip`` is up-to-date. From a terminal (preferably with administrative privileges) use the following command::

		python -m pip install --upgrade pip

Next, DeerLab and all its dependencies can be installed via::

		python -m pip install deerlab

DeerLab will install the following packages:

	* `matplotlib <https://matplotlib.org/>`_ - A comprehensive library for creating data visualizations with Python
	* `memoization <https://pypi.org/project/memoization/>`_ - A powerful caching library for Python
	* `pytest <https://docs.pytest.org/en/stable/>`_ - A Python testing framework
	* `cvxopt <https://cvxopt.org/index.html>`_ - Free software package for convex optimization
	* `numpy <https://numpy.org/>`_ -  Base N-dimensional array package 
	* `scipy <https://www.scipy.org/>`_ - Fundamental library for scientific computing

The installed numerical packages (numpy, scipy, cvxopt) are linked against different BLAS libraries depending on the OS:

	* Windows: linked against the Intel Matrix Kernel Library (MKL)
	* Linux: linked against OpenBLAS
	* Mac: linked against BLAS/LAPACK from the Accelerate framework

Installing older versions
-------------------------

Any DeerLab >0.10.0 release can be installed via pip using the following command matching the x.x.x to the desired version::

		python -m pip install deerlab==x.x.x


Older DeerLab 0.9.x versions, written for MATLAB are available from an `archived repository <https://github.com/JeschkeLab/DeerLab-Matlab>`_. Download and installation instruction for the MATLAB environment are provided in the released documentation. MATLAB releases have been deprecated and no further support is provided.

Installing from source
----------------------

To install DeerLab from the source, first it must be downloaded or cloned from the `DeerLab repository <https://github.com/JeschkeLab/DeerLab>`_. DeerLab and its dependencies can be installed by running the following command on a terminal window to install DeerLab as a static Python package::

		python setup.py install


For developers, in order to install DeerLab but be able to frequently edit the code without having to re-install the package every time use the command::

		python setup.py develop


any further changes made to the source code will then immediate effect.