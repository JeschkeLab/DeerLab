@echo off

set versions=3.6 3.7 3.8 3.9
set platforms=osx-64 linux-32 linux-64 win-32 win-64

:: Get the path to the Anaconda executable
for %%i in (_conda.exe) do (
    set anaconda=%%~p$PATH:i
)
set conda_build=C:%anaconda%conda-bld

echo Activating Anaconda environment...
Powershell.exe -executionpolicy remotesigned C:%anaconda%shell\condabin\conda-hook.ps1 ; conda activate 'C:%anaconda%'

echo Activating Anaconda client...
call anaconda login

echo Building conda package...
:: Delete existing tarball files 
for %%f in (%conda_build%\win-64\*.tar.bz2) do (
    del %%f
)

:: Build  the conda packages for the supported Python versions
for %%v in (%versions%) do (
    call conda-build --python %%v .
)

:: Convert packages to all supported platforms
for %%f in (%conda_build%\win-64\*.tar.bz2) do (
    echo "Converting %%f"
    for %%p in (%platforms%) do (
        call conda-convert --platform %%p %%f  -o %conda_build%
    )
)

:: Upload packages to Anaconda
for %%p in (%platforms%) do (
    for %%f in (%conda_build%\%%p\*.tar.bz2) do (
        call anaconda upload --user JeschkeLab %%f
    )
)

echo "Building conda package finished!"