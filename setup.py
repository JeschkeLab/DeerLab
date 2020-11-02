from setuptools import setup
from setuptools.command.install import install
from setuptools.command.develop import develop
import platform
import subprocess

# Custom OS-specific installation script
def install_dependencies():
#-----------------------------------------------------------------------
    subprocess.run(['pip','install','memoization'],check=True)
    subprocess.run(['pip','install','matplotlib'],check=True)

    # Dependencies with OS-specific BLAS
    if platform.system() == 'Windows':  
        # Install Numpy,SciPy, CVXopt linked to MKL
        subprocess.run(['pip','install','pipwin'],check=True)
        subprocess.run(['pipwin','install','numpy'],check=False)
        subprocess.run(['pipwin','install','scipy'],check=False)
        subprocess.run(['pipwin','install','cvxopt'],check=False)
    else:
        # Install Numpy,SciPy, CVXopt linked to OpenBLAS
        subprocess.run(['pip','install','numpy'],check=True)
        subprocess.run(['pip','install','scipy'],check=True)
        subprocess.run(['pip','install','cvxopt'],check=True)

    subprocess.run(['pip','install','tqdm'],check=True)
    subprocess.run(['pip','install','joblib'],check=True)
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
        install_dependencies()
        subprocess.run(['pip','install','pytest'],check=True) # Install only on development version
#-----------------------------------------------------------------------


setup(
   name='DeerLab',
   version=open('VERSION').read().splitlines()[0],
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
   keywords='data analysis EPR spectroscopy DEER PELDOR',
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
      'Topic :: Scientific/Engineering',
   ]
)