.. _functions:


Functions
---------------------

This is the official documentation for the DeerLab toolbox functions. The following list contains the names of the different function and a brief description of their functionality. The parametric model functions are listed in separate sections.



Modelling
=========================================

This class of functions allows simulation of dipolar signals and their modelling for fitting experimental data.

.. toctree::
    :hidden:
    :glob:
    :maxdepth: 0

    ./api/dipolarkernel
    ./api/aptkernel
    ./api/dipolarsignal
    ./api/whitenoise
    ./api/paramodel
    ./api/mixmodels

.. rst-class:: func-list

=============================================  ============================================================
Function                                         Description
=============================================  ============================================================
:ref:`dipolarkernel`                           Dipolar kernel contructor
:ref:`aptkernel`                               APT kernel contructor
:ref:`dipolarsignal`                           Dipolar signal simulator
:ref:`whitegaussnoise`                         Gaussian white noise generator
:ref:`paramodel`                               Parametric model builder
:ref:`mixmodels`                               Parametric model mixer
=============================================  ============================================================


Analysis
=========================================

This class of functions can be used and/or combined to create fitting routines of dipolar data.

.. toctree::
    :hidden:
    :glob:
    :maxdepth: 0

    ./api/apt
    ./api/backgroundstart
    ./api/fitbackground
    ./api/fitmultigauss
    ./api/fitparamodel
    ./api/fitregmodel
    ./api/obir
    ./api/regoperator
    ./api/sensitivan


.. rst-class:: func-list

=============================================  ============================================================
Function                                         Description
=============================================  ============================================================
:ref:`backgroundstart`                           Background fit start point optimizer
:ref:`fitbackground`                             Background fitting engine
:ref:`fitmultigauss`                             Multi-Gauss fitting engine
:ref:`fitparamodel`                              Parametric model fitting engine
:ref:`fitregmodel`                               Regularization fitting engine
:ref:`obir`                                      Osher-Bregman iterative regularization
:ref:`regoperator`                               Regularization operator constructor
:ref:`apt`                                       Approximate Pake transformation
:ref:`sensitivan`                                Sensitivity analysis engine
=============================================  ============================================================


Pre-Processing
=========================================

This class of functions provide tools for preparing experimental data for analysis. 

.. toctree::
    :hidden:
    :glob:
    :maxdepth: 0

    ./api/correctphase
    ./api/correctzerotime
    ./api/correctscale
    ./api/suppressghost
    ./api/longpass
    ./api/winlowpass

.. rst-class:: func-list

=============================================  ============================================================
Function                                         Description
=============================================  ============================================================
:ref:`correctphase`                              IQ Phase correction
:ref:`correctzerotime`                           Dipolar zero-time correction
:ref:`correctscale`                              Dipolar signal amplitude rescaling
:ref:`suppressghost`                             Ghost-distance suppression
:ref:`longpass`                                  Longpass filtering
:ref:`winlowpass`                                Windowed-lowpass filtering
=============================================  ============================================================



Model Selection
=========================================

This class of functions helps to find an optimal choice of model or model parameters.

.. toctree::
    :hidden:
    :glob:
    :maxdepth: 0

    ./api/selectmodel
    ./api/selregparam
    ./api/regparamrange


.. rst-class:: func-list

=============================================  ============================================================
Function                                         Description
=============================================  ============================================================
:ref:`selectmodel`                              Parametric model selector
:ref:`selregparam`                              Regularization parameter selector
:ref:`regparamrange`                            Regularization parameter range selector
=============================================  ============================================================


Utilities
=========================================

This class of functions provides several tools for quick commands typically required in data processing.


.. toctree::
    :hidden:
    :glob:
    :maxdepth: 0

    ./api/deerload
    ./api/time2freq
    ./api/time2dist
    ./api/noiselevel
    ./api/fftspec
    ./api/prepvalidation


.. rst-class:: func-list

=============================================  ============================================================
Function                                         Description
=============================================  ============================================================
:ref:`deerload`                                  Spectrometer data loader
:ref:`time2freq`                                 Time to frequency axis convertor
:ref:`time2dist`                                 Time to distance axis convertor
:ref:`noiselevel`                                Noise level estimator
:ref:`fftspec`                                   Fast-Fourier trasnform spectrum
:ref:`prepvalidation`                            Full-factorial analysis preparation
=============================================  ============================================================

