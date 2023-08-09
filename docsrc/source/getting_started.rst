.. _beginners_guide:

Getting Started
============================================================

Importing all packages
----------------------

Importing DeerLab
*****************

DeerLab is a Python package. In order to use it, you need to import it. For this, use the import statement: ::

    import deerlab as dl

This makes DeerLab functions accessible via the abbreviated name ``dl``. For example, the function ``deerload`` can be called via ``dl.deerload``. We recommend to use ``dl`` as the standard import abbreviation for DeerLab.

Importing other packages
*************************

Other packages need to be imported as well. The most important one is ::

   import numpy as np  # NumPy: vectors, matrices, linear algebra
   
`NumPy <https://numpy.org/doc/stable/index.html>`_ is the fundamental package for scientific computing in Python. It is a Python library that provides multidimensional arrays (vectors, matrices, etc) and many array functions, including mathematical, logical, shape manipulation, sorting, selecting, I/O, discrete Fourier transforms, basic linear algebra, basic statistical operations, random number generators, and much more.

Most mathematical operations in DeerLab are based on Numpy, and all numerical outputs returned by DeerLab functions are Numpy data types. It is recommendable, to invest a short amount of time to familiarize yourself with some `basic Numpy concepts <https://numpy.org/doc/stable/user/basics.html>`_.

If you have experience with MATLAB, have a look at `Numpy for MATLAB users <https://numpy.org/doc/stable/user/numpy-for-matlab-users.html>`_.

Another important package is `Matplotlib <https://matplotlib.org/>`_, a library that provides plotting capabilities. It contains many modules, of which ``pyplot`` is the most important for basic plotting. Import it with ::

   import matplotlib.pyplot as plt  # Matplotlib: plotting


Python lists and NumPy arrays
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In DeerLab, many functions accept Python lists, NumPy arrays, or both as inputs. While a Python list can contain different data types within a single list, all of the elements in a NumPy array (a so called ndarray) should share the same data type. For example: ::

    a = [1,2,3] # is a list-type
    b = np.array([1,2,3]) # is an ndarray-type

and the elements on such variables can be accessed by their indices in the exact same way: ::

    print(a[0]) # print the first element of the list
    print(b[2]) # print the third element of the ndarray

Note that Python is a 0-indexed language, meaning that the first element of any list, array,... is indexed with 0. 

--------

Loading data and pre-processing
-------------------------------

Loading spectrometer files
***************************

DeerLab provides the function ``deerload`` that can load dipolar EPR data from most spectrometer file formats. It works for 1D and 2D datasets, both real- or complex-valued.

First, determine the location of the spectrometer files you want to load and of the script or Jupyter notebook you are writing is. Let's assume that the script is called ``myscript.py`` and that the data is stored in ``DEERexperiment.DTA`` in the following folder structure: ::

    /home
     |-----experiments
     |      |
     |      |---DEERexperiment.DSC
     |      +---DEERexperiment.DTA
     |
     +-----scripts
            |
            +---myscript.py

From the location of the script, you have two ways to access the data files: using the absolute path, or using the relative path: ::

    filepath = '/home/experiments/DEERexperiment.DTA'   # absolute path
    filepath = '../../experiments/DEERexperiment.DTA'   # relative path

Call ``deerload`` with either of these two paths: ::

    t,Vexp = dl.deerload(filepath)   # load experimental data

The function returns two outputs: the first is the dipolar time-axis of your experiment (a vector of pulse increments), and the second is the raw experimental data as saved by your spectrometer. Here, we store them in variables named ``t`` and ``Vexp``.

Both ``t`` and ``Vexp`` are 1D Numpy arrays with ``N`` elements. To load an additional file, load it into different variables: ::

    filepath1 = '/home/experiments/DEER4p_experiment.DTA'   # absolute path to 1st file
    filepath2 = '/home/experiments/DEER5p_experiment.DTA'   # absolute path to 2nd file
    t1,Vexp2 = dl.deerload(filepath1)   # load 1st set of experimental data
    t2,Vexp2 = dl.deerload(filepath2)   # load 2nd set of experimental data

``deerload`` attempts to return the experiment time-axis ``t`` in units of microseconds, but might not be able to do so for all file formats. For more details about ``deerload`` see the :ref:`reference documentation <deerload>`.

Phase-correction
****************

Experimental dipolar signals are most often acquired in quadrature, with the in-phase and the out-of-phase component stored as the real and the imaginary part of a complex-valued signal. If the out-of-phase components are of no relevance, it is recommendable to perform a phase correction which minimizes the imaginary component and maximizes the real component. If the signal is not complex-valued or the out-of-phase component is important, skip this step. The phase correction function ``correctphase`` takes the complex-valued signal and returns the real-valued phase-corrected dipolar signal: ::

    Vexp = dl.correctphase(Vexp)    # phase correction of experimental data

The correction is based on an optimization approach. This works well in most cases. Should it fail for a specific case, the phase adjustment can also be done manually: ::

    Vexp = np.real(Vexp*np.exp(-1j*phase))    # manual phase correction

