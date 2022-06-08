from sys import executable
import platform
import subprocess

# Custom OS-specific installation script
def install_MKL_dependencies():
#-----------------------------------------------------------------------

    if platform.system() != 'Windows':
        raise SystemError('MKL dependencies are only available on Windows systems.')

    subprocess.run([executable,'-m','pip','uninstall','numpy'],check=False)
    subprocess.run([executable,'-m','pip','uninstall','scipy'],check=False)
    subprocess.run([executable,'-m','pip','uninstall','cvxopt'],check=False)

    # At the moment the dependencies must be handled manually this way
    subprocess.run([executable,'-m','pip','install','git+https://github.com/luisfabib/pipwin'],check=False)

    # Install Numpy,SciPy, CVXopt linked to MKL from Gohlken's repository
    subprocess.run([executable,'-m','pipwin','install','numpy','--filter=mkl','--user'],check=False)
    subprocess.run([executable,'-m','pipwin','install','scipy','--user'],check=False)
    subprocess.run([executable,'-m','pipwin','install','cvxopt','--user'],check=False)
#-----------------------------------------------------------------------    


if __name__ == "__main__":
   print("Installing MKL-linked packages from the Gohlke repository...")
   install_MKL_dependencies()
   print("Finished installing MKL-linked packages from the Gohlke repository.")