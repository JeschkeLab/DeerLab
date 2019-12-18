## Compiling the documentation

The DeerAnalysis documentation/website is written in restructred-text format and built with Sphinx. In order to compile the source files of the documentation, several packages/programs are required. This is a summary of the steps to be taken before being able to compile the documentation.

### Requirements: 
  * Pyhton3
  * Sphinx 1.8 (or older)
  * Matlab Sphinx-domain
  * Read-the-Docs Sphinx-theme
	
### Installation:

In order to compile the documentation the following steps must be followed:

1) Install and setup python environment from https://www.python.org/

2) Install sphinx

    * From python console
            
                pip install sphinx<=1.8.0
    * From DOS console
    
                python -m pip install sphinx<=1.8.0

3) Install Read-the-Docs Sphinx theme

    * From python console
    
            pip install sphinx_rtd_theme
    * From DOS console
    
            python -m pip install sphinx_rtd_theme

4) Install Matlab sphinx-domain

    * From python console
    
            pip install sphinxcontrib-matlabdomain
    * From DOS console
    
            python -m pip install sphinxcontrib-matlabdomain
	
5) Install HTTP sphinx-domain

    * From python console
    
            pip install sphinxcontrib-httpdomain
    * From DOS console		
    
            python -m pip install sphinxcontrib-httpdomain
		
6) Download and install dvissvgm from https://dvisvgm.de/Downloads/
		
7) From /DeerAnalysis/docsrc/ run the batch script

	    #To compile using the cached data
            make    
	    #To compile from scratch
            make clean
