from setuptools import setup

setup(
    name='DeerLab',
    version=open('VERSION').read().splitlines()[0].replace("v", ""),
    author='Luis Fábregas Ibáñez , Stefan Stoll and other contributors',
    package_dir={'deerlab': 'deerlab'},
    packages=['deerlab','deerlab'],
    url='https://github.com/JeschkeLab/DeerLab',
    project_urls = {
        'Documentation': 'https://jeschkelab.github.io/DeerLab/',
        'Source': 'https://github.com/JeschkeLab/DeerLab',
    },
    python_requires='>=3.9',
    license='LICENSE.txt',
    include_package_data = True,
    keywords='data analysis modeling least-squares EPR spectroscopy DEER PELDOR'.split(),
    description='Comprehensive package for data analysis of dipolar EPR spectroscopy',
    long_description=open('README.md', encoding='utf-8').read(),
    long_description_content_type="text/markdown",
    install_requires = [
                        'numpy>=2.0',
                        'cvxopt>=1.0.0',
                        'scipy>=1.11',
                        'joblib>=1.0.0',
                        'dill>=0.3.0',
                        'tqdm>=4.51.0',
                        'matplotlib>=3.6.0',
                        'memoization>=0.3.1',
                        'pytest>=6.2.2',
                        'setuptools>=53.0.0',
                        'numexpr>=2.7.3',
                        'quadprog>=0.1.11; python_version <= "3.10"',
                        ],
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Operating System :: Microsoft :: Windows',
        'Operating System :: MacOS',
        'Operating System :: POSIX :: Linux',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'Programming Language :: Python :: 3.12',
        'Topic :: Scientific/Engineering :: Chemistry',
    ]
)
