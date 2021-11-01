.. _beginners_guide:

Getting Started
============================================================

This is the introductory guide to DeerLab for any dipolar electron paramagnetic resonance (EPR) spectroscopy applications.

--------

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

    t,V = dl.deerload(filepath)   # load experimental data

The function returns two outputs: the first is the dipolar time-axis of your experiment (a vector of pulse increments), and the second is the raw experimental data as saved by your spectrometer. Here, we store them in variables named ``t`` and ``V``.

Both ``t`` and ``V`` are 1D Numpy arrays with ``N`` elements. To load an additional file, load it into different variables: ::

    filepath1 = '/home/experiments/DEER4p_experiment.DTA'   # absolute path to 1st file
    filepath2 = '/home/experiments/DEER5p_experiment.DTA'   # absolute path to 2nd file
    t1,V2 = dl.deerload(filepath1)   # load 1st set of experimental data
    t2,V2 = dl.deerload(filepath2)   # load 2nd set of experimental data

``deerload`` attempts to return the experiment time-axis ``t`` in units of microseconds, but might not be able to do so for all file formats. For more details about ``deerload`` see the :ref:`reference documentation <deerload>`.

Phase-correction
****************

Experimental dipolar signals are most often aquired in quadrature, with the in-phase and the out-of-phase component stored as the real and the imaginary part of a complex-valued signal. If the out-of-phase components are of no relevance, it is recommendable to perform a phase correction which minimizes the imaginary component and maximizes the real component. If the signal is not complex-valued or the out-of-phase component is important, skip this step. The phase correction function ``correctphase`` takes the complex-valued signal and returns the real-valued phase-corrected dipolar signal: ::

    V = dl.correctphase(V)    # phase correction of experimental data

The correction is based on an optimization approach. This works well in most cases. Should it fail for a specific case, the phase adjustment can also be done manually: ::

    V = np.real(V*np.exp(-1j*phase))    # manual phase correction

---------------

Dipolar modelling
-------------------------

DeerLab provides a very flexible framework to model dipolar signals originating from any dipolar EPR spectroscopy experiments. Choosing a model that properly describes your sample and experiment is of paramount importance. The DeerLab function ``dipolarmodel`` already defines the core model structure based on dipolar pathways, with the following components to be chosen:     

* **Distance range**: Also called the interspin distance axis, is the range of distances where the distribution is defined. 

* **Distribution model**: Describes the intra-molecular distance distribution in either a parametric (e.g. a Gaussian distribution) or a non-parametric way. 

* **Background model**: Describes the dipolar background signal arising from the inter-molecular contributions. 

* **Number of pathways**: Sets the number of dipolar pathways contributing to the dipolar signal.

For each of these four components, a choice needs to be made: 

Choosing a distance range
*************************

The distance range :math:`[r_\mathrm{min},r_\mathrm{max}]` is an important choice, as any distance distribution is truncated to this range, i.e. :math:`P(r)=0` for :math:`r<r_\mathrm{min}` and :math:`r>r_\mathrm{max}`. The lower limit of the distance range is determined by the bandwidth of the pulses, and also by the time increment. Typically, 1.5 nm is a reasonable choice. The upper limit depends on the distances in your sample. The number of points in ``r`` is usually set to a certain resolution (typically 0.01-0.05nm). Such a distance-axis is usually defined as ``r`` is most easily defined using the ``linspace`` function from NumPy: ::

    r = np.linspace(1.5,6.5,100)  # define distance range from 1.5nm to 6.5nm with a resolution of 0.05nm

Choosing a distribution model
******************************

A non-parametric distribution is specified by setting the choice of ``Pmodel`` keyword in ``dipolarmodel`` to ``None``. In a non-parametric distribution, each element :math:`P_i` of the distribution is a linear parameter. Non-parametric distributions are obtained via methods such as Tikhonov regularization. If there are reasons to believe that the distance distribution has a specific shape (e.g. Gaussian, Rice, random-coil, etc.), or if there is very little information in the data, use a parametric distance distribution model from the :ref:`list of available models<modelsref_dd>`.

Choose a background model
*************************

Typically, a background model of a homogenous 3D distribution of spins is appropriate. The associated parametric model function is :ref:`bg_hom3d`. In some cases, depending on the properties of your sample, other background models might be needed, such as backgrounds arising from distributions of spins in fractal dimensions or when  accounting for volume-exclusion effects. In such cases, use the associated parametric background models from the :ref:`list of available models<modelsref_bg>`. If there is no inter-molecular background in your sample, or it is negligible, set the background model to ``None``.



Choosing the number of dipolar pathways
*************************************** 

This decision should be based on the experiment you used to acquire the data and the type of pulses you used. In the case of 4-pulse DEER data, when analyzing a standard 4-pulse DEER signal without 2+1 component at the end a single pathway suffices. If the 2+1 component (appearing at the right edge of the time trace) is present, then it should be fitted as well, including its counterpart appearing at negative times, making a total of three dipolar pathways. Experiments such as 5-pulse DEER typically require at least two dipolar pathways to be properly modelled. 


Using experimental pulse delays
******************************** 

The dipolar pathways of a newly constructed dipolar model are initialized at arbitrary refocusing times and fully unconstrained. The refocusing times can be strongly constrained by knowing the experimental pulse sequence delays used to acquire the data. If the experiment used to acquire the data is known, as well as its pulse delays, then it is strongly recommended do so.  
 
DeerLab provides a selection of experimental information generators for some of the most widely employed experimental methods (see the of :ref:`list of available experiments <modelsref_ex>`). These are functions that take the pulse sequence delays, and return an ``ExperimentInfo`` object. This can be passed to the ``dipolarmodel`` function via the ``experiment`` keyword argument, to incorporate the experiment information on the model and constrain some of its parameters. 

Constructing the dipolar model 
*******************************

Once all the decisions above have been made, the dipolar model can be constructed using the ``dipolarmodel`` function. The models that have an associated parametric function, e.g. ``bg_hom3d``, must be passed directly as inputs to ``dipolarmodel``. In Python, functions can be passed as inputs to other functions.  See the :ref:`details <dipolarmodel>` on ``dipolarmodel`` for more information. 

Example: Single-pathway 4-pulse DEER model
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
For example, a 4pDEER signal with non-parametric distance distribution and homogenous 3D background can be constructed using ::

    expinfo = dl.ex_4pdeer(tau1=0.5, tau2=5.5)
    Vmodel = dl.dipolarmodel(t, r, Pmodel=None, Bmodel=dl.bg_hom3d, npathways=1, experiment=expinfo) 

By default, the function ``dipolarmodel`` assumes a non-parametric distance distribution, a homogenous 3D background and a single pathway. Thus the above is equivalent to ::

    expinfo = dl.ex_4pdeer(tau1=0.5, tau2=5.5)
    Vmodel = dl.dipolarmodel(t, r, experiment=expinfo) 




Example: Two-pathway 5-pulse DEER model
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
For example, a 5pDEER signal with non-parametric distance distribution and homogenous 3D background can be constructed using ::

    expinfo = dl.ex_5pdeer(tau1=0.5, tau2=5.5, tau3=0.2)
    Vmodel = dl.dipolarmodel(t, r, Pmodel=None, Bmodel=dl.bg_hom3d, npathways=2, , experiment=expinfo)

Manipulating the model
***********************

A full summary of the constructed model(s) can be inspected by printing the model object ::

    >>> print(Vmodel)
    Model information 
    -----------------

    Model description: Dipolar signal model
    Model call signature: (mod,reftime,conc,P)
    Constants: []

    Parameter Table 
    ---------------

    ============ ========= ========== =========== ======== ========== ==========================
        Name       Lower     Upper       Type      Frozen    Units      Description  
    ============ ========= ========== =========== ======== ========== ==========================
      mod           0         1          nonlin      No                 Modulation depth
      reftime       -inf      inf        nonlin      No       μs        Refocusing time
      conc          0.01      5e+03      nonlin      No       μM        Spin concentration
      P             0         inf        linear      No       None      Non-parametric distance distribution
    ============ ========= ========== =========== ======== ========== ==========================


From this point on, the model can be modified, manipulated and expanded freely as any other DeerLab model. Check out the :ref:`modelling guide <modelling_guide>` for more details and instructions on model manipulation.

Fitting
-------
Next, the model ``Vmodel`` can be fitted to the experimental data ``V`` by calling the ``fit`` function: ::

    result = dl.fit(Vmodel,V)  # Fit the model to the experimental data


After ``fit`` has found a solution, it returns an object that we assigned to ``result``. This object contains fields with all quantities of interest with the fit results, such as the fitted model and parameters, goodness-of-fit statistics, and uncertainty information. Check out the :ref:`fitting guide <fitting_fitresult>` for more details on the quantities provided in ``result``.


Displaying the results
**********************

For just a quick display of the results, you can use the ``plot()`` method of the ``fit`` object that will display a figure with you experimental data, the corresponding fit including confidence bands. :: 

    fitresults.plot() # display results


.. image:: ./images/beginners_guide1.png
   :width: 450px

These confidence bands are covariance-based and might represent an overestimation of the true uncertainty on the results (see :ref:`uncertainty <uncertainty>` for further details). It is important to always report confidence bands with fitted distance distributions.

The ``fitresults`` output contains additional information, for example:

    * ``fit.V``, ``fit.B``, and ``fit.P`` contain the arrays of the fitted dipolar signal, background, and distance distribution, respectively. 
    * ``fit.exparam``, ``fit.bgparam``, and ``fit.ddparam`` contain the arrays of fitted model parameters for the experiment, background, and distribution models. 
    * ``fit.scale`` contains the fitted overall scale of the dipolar signal.

In addition to the distance distribution fit, it is important to check and report the fitted model parameters and their uncertainties. While this can be computed manually, a summary can be easily requested by enabling the ``verbose`` option of ``fitmodel``. By using ::

    fit = dl.fitmodel(V,t,r,'P',dl.bg_hom3d,dl.ex_4pdeer,verbose=True)  # 4pDEER fit and report parameter fits

after the function has fitted your data, it will print a summary the results, including goodness-of-fit estimators
and fitted parameters with uncertainties. Here is an example output

.. code-block:: text

    -----------------------------------------------------------------------------------------
    Goodness of fit
    Vexp[0]: 𝛘2 = 25.510184  RMSD  = 1.953580e+07
    -----------------------------------------------------------------------------------------
    Fitted parameters and 95%-confidence intervals
    Vfit[0]:
    V0:  3.551e+07  Signal scale (arb.u.)
    bgparam[0]:   145.3121342  (111.0809911, 179.5432773)  Concentration of pumped spins (μM)
    exparam[0]:   0.4066627  (0.3630338, 0.4502916)  Modulation depth ()
    -----------------------------------------------------------------------------------------

where there are no distribution parameters (``ddparam``) due to the distribution model being non-parametric. 


Exporting the figure and the data
*********************************

After completing the fit, you might want to export the figure with the fit. Here is one way to do it: ::

    figure = fit.plot()                       # get figure object
    figure.savefig('DEERFig.png', dpi=600)    # save figure as png file
    figure.savefig('DEERFig.pdf')             # save figure as pdf file

To export the fitted distance distribution for plotting with another software, save it in a simple text file ::

    np.savetxt('distancedistribution.txt', np.asarray((r, fit.P, *fit.Puncert.ci(95).T)).T)

The generated file contain four columns: the distance axis, the distance distributions, and the upper and lower confidence bounds. The ``.T`` indicate array transposes, which are used to get the confidence bands into the column format for saving.

To export the fitted time-domain trace, use similarly ::

    np.savetxt('timetrace.txt', np.asarray((t, V, fit.V, *fit.Vuncert.ci(95).T)).T)

------------

Summary
--------

Here is an example script to load experimental time trace, pre-process it, and fit a 4-pulse DEER model with a non-parametric distance distribution:  ::

    import numpy as np
    import deerlab as dl

    # Optional, if experimental delays known
    expinfo = dl.ex_4pdeer(tau1=0.5,tau2=5.5)

    # Data import and pre-processing
    filepath = '/home/experiments/DEERexperiment.DTA' # File path
    t,Vexp = dl.deerload(filepath) # Load experimental data
    Vexp = dl.correctphase(Vexp) # Phase correction 

    # Distance range
    r = np.linspace(1.5,6.5,100) # Define distance range from 1.5nm to 6nm with a resolution of 0.05nm
    
    # Construct the dipolar model 
    Vmodel = dl.dipolarmodel(t,r) # Non-parametric P(r), homogenous 3D background, single-pathway

    # Fit the model to the data
    result = dl.fit(Vmodel,Vexp)

    # Print figure
    figure = result.plot()
    figure.savefig('DEERfig.pdf')
