#  Global environment: Python 3.6.0 (64-bit, WinOS) 

# Virtual environment preparation 
# Launches the virtual environment for DeerLab and prepares the activation scripts to link the MKL binaries
import fileinput
import subprocess
import venv

# Create virtual environment
subprocess.run(["python","-m", "venv",".venv"])

# Edit the activation scripts for the virtual environment to ensure that 
# the Intel MKL binaries are included in the search path

# Batch script
for line in fileinput.input('./.venv/Scripts/Activate.ps1', inplace=True):
    if '$env:PATH = "$env:VIRTUAL_ENV\\Scripts;$env:PATH"' in line:
        line = line.replace(line,'$env:PATH = "$env:VIRTUAL_ENV\\Scripts;$env:VIRTUAL_ENV\\Library;$env:VIRTUAL_ENV\\Library\\bin;$env:PATH"\n')
    print(line)
# Powershell script
for line in fileinput.input('./.venv/Scripts/activate.bat', inplace=True):
    if 'set "PATH=%VIRTUAL_ENV%\\Scripts;' in line:
        line = line.replace(line,'set "PATH=%VIRTUAL_ENV%\\Scripts;%VIRTUAL_ENV%\\Library;%VIRTUAL_ENV%\\Library\\bin;%PATH%"\n')
    print(line)



