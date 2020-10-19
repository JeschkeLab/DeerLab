## Compiling the documentation

The DeerLab documentation/website is written in restructured-text (RST) format and built using Sphinx. In order to compile the source files of the documentation, several packages/programs are required. This is a summary of the steps to be taken before being able to compile the documentation.

### Requirements: 
  * Pyhton3
  * Sphinx 1.8 (or older)
  * numpydoc
  * Sphinx-Gallery
  * Read-the-Docs Sphinx-theme
  * DeerLab
  * dvissvgm

### Installation:

In order to compile the documentation the following steps must be followed:

1) Install and setup python environment from https://www.python.org/

2) Install DeerLab (see installation instructions)

2) Install Sphinx

        pip install sphinx==1.8.0

3) Install Read-the-Docs Sphinx theme
    
        pip install sphinx_rtd_theme

5) Install HTTP sphinx-domain
    
        pip install sphinxcontrib-httpdomain

6) Install numpydoc
    
        pip install numpydoc

7) Install Sphinx-Gallery
    
        pip install sphinx-gallery
        
8) Install M2R-2
    
        pip install m2r2

9) Download and install dvissvgm from https://dvisvgm.de/Downloads/
		
### Runnning the Sphinx builder


To build the documentation from the source, call the Makefile or make.bat scripts:

        #To compile using the cached data
        ./docsrc/make.bat    
        #To compile from scratch
        ./docsrc/make.bat clean    
