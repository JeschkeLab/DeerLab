COMPILATION OF DEERANALYSIS2 DOCUMENTATION

Requirements: 
	Pyhton3
	Sphinx 1.8 (or older)
	Matlab Sphinx-domain
	Read-the-Docs Sphinx-theme
	
Installation:

0) Install and setup python environment from

    https://www.python.org/

1) Install sphinx

    - From python console
            
                pip install sphinx<=1.8.0
    - From DOS console
    
                python -m pip install sphinx<=1.8.0

3) Install Read-the-Docs Sphinx theme

    - From python console
    
            pip install sphinx_rtd_theme
    - From DOS console
    
            python -m pip install sphinx_rtd_theme

2) Install Matlab sphinx-domain

    - From python console
    
            pip install sphinxcontrib-matlabdomain
    - From DOS console
    
            python -m pip install sphinxcontrib-matlabdomain
	
3) Install HTTP sphinx-domain
	
    -From python console
    
        pip install sphinxcontrib-httpdomain
    -From DOS console		
    
        python -m pip install sphinxcontrib-httpdomain
		
4) Download and install dvissvgm from

	https://dvisvgm.de/Downloads/
		
5) From /DeerAnalysis/docsrc/ run the batch script

        make html