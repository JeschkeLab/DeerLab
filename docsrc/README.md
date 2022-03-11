## Compiling the documentation

The DeerLab documentation/website is written in restructured-text (RST) format and built using [Sphinx](https://www.sphinx-doc.org/). In order to compile the source files of the documentation, several packages and programs are required. This is a summary of the steps to be taken before being able to compile the documentation.

### Installating required packages:

In order to compile the documentation the following steps must be followed:

1) Install and setup [Python](https://www.python.org/)

2) Install DeerLab from source (see [installation instructions](https://jeschkelab.github.io/DeerLab/installation.html))

3) Install Sphinx and Sphinx extensions

        pip install sphinx    
        pip install pydata-sphinx-theme
        pip install sphinxcontrib-httpdomain
        pip install sphinx-gallery
        pip install sphinxcontrib-ghcontributors
        pip install sphinx-issues
        pip install sphinx-copybutton 
        pip install numpydoc

5) Download and install [dvisvgm](https://dvisvgm.de/Downloads/).
		
### Runnning the Sphinx builder

To build the documentation from the source, call the `Makefile` or `make.bat` scripts:

        # To compile using the cached data
        ./docsrc/make.bat
        # To compile from scratch
        ./docsrc/make.bat clean    

If the `sphinx-build` command is not found, the documentation can be built with the following command 

        python -m sphinx.cmd.build -d ./cache -b html ./source ../docs

To skip running the examples with the Sphinx Gallery extension (can take very long on a first build), use the option `-D  plot_gallery=0`, e.g.

        python -m sphinx.cmd.build -E -b html ./source ../docs -D  plot_gallery=0
