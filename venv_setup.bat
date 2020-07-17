:: ===========================================
::  DEERLAB - VIRTUAL ENVIRONMENT CONSTRUCTOR 
:: ===========================================
::  Prepares all packages and dependencies, 
::  including Intel MKL libraries linked to 
::  Scipy and Numpy. 
:: ===========================================
:: REQUIRES: Python 3.6
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
call python -m pip install --upgrade setuptools
:: Install all packages in the requirements list
for /f "tokens=*" %%a in (requirements.txt) do (
    IF "%%a"=="intel-scipy" (
        :: The intel-scipy distribution includes the intel-numpy package
        call python -m pip install %%a
     ) ELSE (
        :: Install all other dependencies without updating or reinstalling numpy
        call python -m pip install --upgrade-strategy only-if-needed %%a
     ) 
)
:: Setup DeerLab
call python setup.py develop
