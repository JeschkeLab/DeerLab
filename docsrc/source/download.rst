Installation
======================

Get the latest version of DeerAnalysis from the `GitHub repository <https://github.com/luisfabib/DeerAnalysis2>`_.


Requirements
---------------
The application programming interface (API) of DeerAnalysis requires the following products:

    *  `MATLAB <https://ch.mathworks.com/products/matlab.html>`_
    
    *  `Optimization Toolbox <https://ch.mathworks.com/products/optimization.html>`_

Setup
---------------
In order for MATLAB to access the DeerAnalysis API functions, the path to DeerAnalysis installation folder must be set.

Method A) Add DeerAnalysis path via MATLAB's IDE
	1) On the ``Home`` tab, in the ``Environment`` section, click ``Set Path``. 
	2) Click ``Add with Subfolders...`` and select the ``DeerAnalysis\functions`` directory. 
	3) Click ``Save`` to save the current MATLAB search path and exit via ``Close``.


Method B) Adding DeerAnalysis path at startup
	1) Open (or create) the ``startup.m`` file in the default ``\MATLAB`` directory.
	2) Add the following lines of code:


	 .. code-block:: matlab

		 addpath('mypath/DeerAnalysis/functions')
		 addpath('mypath/DeerAnalysis/functions/models')
		 
	3) Save ``startup.m`` and restart MATLAB.
