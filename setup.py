from setuptools import setup

setup(
   name='DeerLab',
   version='0.10.0',
   author='Luis Fabregas, Stefan Stoll and other contributors',
   packages=['deerlab', 'deerlab.test'],
   url='https://jeschkelab.github.io/DeerLab/',
   license='LICENSE.txt',
   description='An awesome package that does something',
   long_description=open('README.md').read(),
   install_requires=[
       "pytest",
       "numpy",
       "pandas",
       "scipy"
   ],
)