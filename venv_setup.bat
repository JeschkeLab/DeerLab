:: ===========================================
::  DEERLAB - VIRTUAL ENVIRONMENT CONSTRUCTOR 
:: ===========================================
::  Prepares all packages and dependencies, 
::  including Intel MKL libraries linked to 
::  Scipy and Numpy. 
:: ===========================================
::  Can be executed via the terminal: 
::      .\venv_setup.bat
::  Afterwards the virtual environment can be
::  activated via the command:
::      ./.venv/Scripts/activate
:: ===========================================

:: Construct and configure virtual environment
call python -m venv_prep
:: Activate the environment and return
call ./.venv/Scripts/activate
:: pip manager is crude, upgrade 
call python -m pip install --upgrade pip
:: Install package dependencies of DeerLab
call python -m pip install -r requirements.txt
:: Due to dependency overlap, some packages might overlap 
:: with intel-numpy, delete all numpy versions
call python -m pip uninstall numpy -y
call python -m pip uninstall numpy -y
:: Re-install the proper Intel numpy version
call python -m pip install intel-numpy
:: Setup DeerLab
call python setup.py develop
