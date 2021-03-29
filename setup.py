from setuptools import setup
from setuptools.command.install import install
from setuptools.command.develop import develop
from sys import executable
import platform
import subprocess

# Custom OS-specific installation script
def install_dependencies(develop_mode=False):
#-----------------------------------------------------------------------

    # At the moment the dependencies must be handled manually this way due to 
    # the dependency of Windows systems on the external Gohlke binaries. This
    # also means that the package must be distributed as a source tarball instead
    # of the commonly used wheels. 
    # As soon as Numpy, Scipy and CVXOpt start distributing their own binaries for 
    # Windows systems this will no longer be necessary and all dependencies will be 
    # able to be handled by pip automatically. 

    # DeerLab dependencies on the PyPI
    dependencies = ["memoization","matplotlib","tqdm","joblib"] 

    if platform.system() == 'Windows':
        # Use own forked pipwin repo for self-patched fixes 
        dependencies += ["git+https://github.com/luisfabib/pipwin"]
    else:
        dependencies += ["numpy","scipy","cvxopt"]
    if develop_mode:
        dependencies += ["pytest"]
        
    # Commands for pip installation of dependencies
    pipinstall = [executable,'-m','pip','install']
    userlevel = ['--user']

    # Easier to Ask Forgiveness than Permission Principle: 
    # try to install dependencies system-wide, if fails then at user-level
    for dependency in dependencies:
        try:
            subprocess.run(pipinstall+[dependency],check=True)
        except:
            subprocess.run(pipinstall+userlevel+[dependency],check=False)

    # Download and install pre-compiled binaries for Windows-systems
    if platform.system() == 'Windows':  
        # Refresh the pipwin cache to get the latest repo status
        subprocess.run(['pipwin','refresh'],check=False)
        # Install Numpy,SciPy, CVXopt linked to MKL from Gohlken's repository
        subprocess.run(['pipwin','install','numpy','--filter=mkl'],check=False)
        subprocess.run(['pipwin','install','scipy'],check=False)
        subprocess.run(['pipwin','install','cvxopt'],check=False)
#-----------------------------------------------------------------------    

class install_routine(install):
#-----------------------------------------------------------------------
    """Customized setuptools install command"""
    def run(self):
        install.run(self)
        install_dependencies()
#-----------------------------------------------------------------------

class develop_routine(develop):
#-----------------------------------------------------------------------
    """Customized setuptools install command"""
    def run(self):
        develop.run(self)
        install_dependencies(develop_mode=True)
#-----------------------------------------------------------------------

setup(
    name='DeerLab',
    version=open('VERSION').read().splitlines()[0].replace("v", ""),
    author='Luis Fábregas Ibáñez , Stefan Stoll and other contributors',
    package_dir={'deerlab': 'deerlab',
                'deerlab.utils': 'deerlab/utils'},
    packages=['deerlab','deerlab.utils'],
    url='https://github.com/JeschkeLab/DeerLab',
    project_urls = {
        'Documentation': 'https://jeschkelab.github.io/DeerLab/',
        'Source': 'https://github.com/JeschkeLab/DeerLab',
    },
    python_requires='>=3.6',
    license='LICENSE.txt',
    include_package_data = True,
    keywords='data analysis EPR spectroscopy DEER PELDOR'.split(),
    description='Comprehensive package for data analysis of dipolar EPR spectroscopy',
    long_description=open('README.md').read(),
    long_description_content_type="text/markdown",
    cmdclass={
        'install': install_routine,
        'develop': develop_routine
    },
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Operating System :: Microsoft :: Windows',
        'Operating System :: MacOS',
        'Operating System :: POSIX :: Linux',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Topic :: Scientific/Engineering',
    ]
)