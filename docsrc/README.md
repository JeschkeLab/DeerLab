## Compiling the documentation

The DeerLab documentation/website is written in restructured-text (RST) format and built using Sphinx. In order to compile the source files of the documentation, several packages and programs are required. This is a summary of the steps to be taken before being able to compile the documentation.

### Installation:

In order to compile the documentation the following steps must be followed:

1) Install and setup Python environment from https://www.python.org/

2) Install DeerLab (see installation instructions)

3) Install Sphinx

        pip install sphinx

4) Install PyData Sphinx theme
    
        pip install pydata-sphinx-theme

5) Install HTTP sphinx-domain
    
        pip install sphinxcontrib-httpdomain

6) Install numpydoc
    
        pip install numpydoc

7) Install Sphinx-Gallery
    
        pip install sphinx-gallery
        
8) Install ghcontributors

        pip install sphinxcontrib-ghcontributors

9) Install Sphinx-Issues
    
        pip install sphinx-issues

10) Install Sphinx-Copy Button

        pip install sphinx-copybutton 

11) Download and install dvissvgm from https://dvisvgm.de/Downloads/
		
### Runnning the Sphinx builder

To build the documentation from the source, call the Makefile or make.bat scripts:

        # To compile using the cached data
        ./docsrc/make.bat    
        # To compile from scratch
        ./docsrc/make.bat clean    

If the `sphinx-build` command is not found, the documentation can be built with the following command 

        python -m sphinx.cmd.build -d ./cache -b html ./source ../docs

To skip running the examples with the Sphinx Gallery extension (can be very long on a first build) use the option `-D  plot_gallery=0`, e.g.

        python -m sphinx.cmd.build -E -b html ./source ../docs -D  plot_gallery=0  