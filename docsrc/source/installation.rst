Installation
======================


Requirements
---------------
The application programming interface (API) of DeerLab requires the following products:

    *  `MATLAB <https://ch.mathworks.com/products/matlab.html>`_


To access additional optional functionality, the following products may be required:

    *  `Optimization Toolbox <https://ch.mathworks.com/products/optimization.html>`_


-----------------------


Setup
---------------
In order for MATLAB to access the DeerLab API functions, the path to DeerLab installation folder must be set.



Method A) Add DeerLab path via MATLAB's IDE
	1) On the ``Home`` tab, in the ``Environment`` section, click ``Set Path``. 
	
		.. image:: ./images/installation1.png
		
	2) Click ``Add with Subfolders...`` and select the ``DeerLab\functions`` directory. 
	
		.. image:: ./images/installation2.png
		
	3) Click ``Save`` to save the current MATLAB search path and exit via ``Close``.


-----------------------

Method B) Add DeerLab path at startup
	1) Open (or create) the ``startup.m`` file in the default ``\MATLAB`` directory.
	2) Add the following lines of code:


	 .. code-block:: matlab

		 addpath('mypath/DeerLab/functions')
		 addpath('mypath/DeerLab/functions/models')
		 
	3) Save ``startup.m`` and restart MATLAB.
