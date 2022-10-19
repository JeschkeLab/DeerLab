from setuptools import setup

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
    python_requires='>=3.8',
    license='LICENSE.txt',
    include_package_data = True,
    keywords='data analysis modelling least-squares EPR spectroscopy DEER PELDOR'.split(),
    description='Comprehensive package for data analysis of dipolar EPR spectroscopy',
    long_description=open('README.md', encoding='utf-8').read(),
    long_description_content_type="text/markdown",
    install_requires = [
                        'numpy<=1.22.4',
                        'cvxopt>=1.0.0',
                        'scipy>=1.6.3',
                        'joblib>-1.0.0',
                        'dill>=0.3.0',
                        'tqdm>=4.51.0',
                        'matplotlib>=3.3.4',
                        'memoization>=0.3.1',
                        'pytest>=6.2.2',
                        'setuptools>=53.0.0',
                        'numexpr >=2.7.3',
                        ],
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Operating System :: Microsoft :: Windows',
        'Operating System :: MacOS',
        'Operating System :: POSIX :: Linux',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Topic :: Scientific/Engineering',
    ]
)