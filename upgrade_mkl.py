from sys import executable
import platform
import subprocess
import os

MKL_LINKED = os.environ.get('MKL_LINKED', False)

# Custom OS-specific installation script
def install_MKL_dependencies():
#-----------------------------------------------------------------------

    if platform.system() != 'Windows':
        raise SystemError('MKL dependencies are only available on Windows systems.')

    # At the moment the dependencies must be handled manually this way
    subprocess.run([executable,'-m','pip','install','git+https://github.com/luisfabib/pipwin'],check=False)

    try:
        # Install Numpy,SciPy, CVXopt linked to MKL from Gohlken's repository
        subprocess.run([executable,'-m','pipwin','install','numpy','--filter=mkl'],check=False)
        subprocess.run([executable,'-m','pipwin','install','scipy'],check=False)
        subprocess.run([executable,'-m','pipwin','install','cvxopt'],check=False)
    except: 
        subprocess.run([executable,'-m','pipwin','install','numpy','--filter=mkl','--user'],check=False)
        subprocess.run([executable,'-m','pipwin','install','scipy','--user'],check=False)
        subprocess.run([executable,'-m','pipwin','install','cvxopt','--user'],check=False)
#-----------------------------------------------------------------------    


if __name__ == "__main__":
   print("Installing MKL-linked packages from the Gohlke repository...")
   install_MKL_dependencies()
   print("Finished installing MKL-linked packages from the Gohlke repository.")