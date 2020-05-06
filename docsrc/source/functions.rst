.. _functions:


Functions
---------------------

This is the official documentation for the DeerLab toolbox functions. The following list contains the names of the different function and a brief description of their functionality. The parametric model functions are listed in separate sections.



Modeling
=========================================

This class of functions allows simulation of dipolar signals and their modeling for fitting experimental data.

.. toctree::
    :hidden:
    :glob:
    :maxdepth: 0

    ./functions/dipolarkernel
    ./functions/aptkernel
    ./functions/dipolarsignal
    ./functions/dipolarbackground
    ./functions/whitenoise
    ./functions/paramodel
    ./functions/mixmodels

.. rst-class:: func-list

=============================================  ============================================================
Function                                         Description
=============================================  ============================================================
:ref:`dipolarkernel`                           Dipolar kernel constructor
:ref:`dipolarbackground`                       Multi-pathway background constructor
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

    ./functions/fitsignal
    ./functions/apt
    ./functions/backgroundstart
    ./functions/fitbackground
    ./functions/fitmultimodel
    ./functions/fitparamodel
    ./functions/fitregmodel
    ./functions/obir
    ./functions/regoperator
    ./functions/bootan
    ./functions/sensitivan


.. rst-class:: func-list

=============================================  ============================================================
Function                                         Description
=============================================  ============================================================
:ref:`fitsignal`                                 Full model signal fitting engine
:ref:`bootan`                                    Bootstrap uncertainty analysis
:ref:`fitmultimodel`                             Multi-component model fitting engine
:ref:`fitparamodel`                              Parametric model fitting engine
:ref:`fitregmodel`                               Regularization fitting engine
:ref:`obir`                                      Osher-Bregman iterative regularization
:ref:`regoperator`                               Regularization operator constructor
=============================================  ============================================================


Pre-Processing
=========================================

This class of functions provide tools for preparing experimental data for analysis. 

.. toctree::
    :hidden:
    :glob:
    :maxdepth: 0

    ./functions/correctphase
    ./functions/correctzerotime
    ./functions/correctscale
    ./functions/longpass
    ./functions/winlowpass

.. rst-class:: func-list

=============================================  ============================================================
Function                                         Description
=============================================  ============================================================
:ref:`correctphase`                              IQ Phase correction
:ref:`correctzerotime`                           Dipolar zero-time correction
:ref:`correctscale`                              Dipolar signal amplitude rescaling
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

    ./functions/selectmodel
    ./functions/selregparam
    ./functions/regparamrange


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

    ./functions/deerload
    ./functions/time2freq
    ./functions/time2dist
    ./functions/noiselevel
    ./functions/fftspec
    ./functions/prepvalidation


.. rst-class:: func-list

=============================================  ============================================================
Function                                         Description
=============================================  ============================================================
:ref:`deerload`                                  Spectrometer data loader
:ref:`time2freq`                                 Time to frequency axis converter
:ref:`noiselevel`                                Noise level estimator
:ref:`fftspec`                                   Fast-Fourier transform spectrum
=============================================  ============================================================

---------------------------

Obsolescent Functions
""""""""""""""""""""""


This class of functions provides tools for reproducing obsolete analysis methods or workflows encountered in older software (e.g. DeerAnalysis). These functions have become obsolete and are not recommended for routine data analysis.

.. rst-class:: func-list

=============================================  ============================================================
Function                                         Description
=============================================  ============================================================
:ref:`time2dist`                                 Time to distance axis converter
:ref:`backgroundstart`                           Background fit start point optimizer
:ref:`fitbackground`                             Background fitting engine
:ref:`aptkernel`                                 APT kernel constructor
:ref:`apt`                                       Approximate Pake transformation
=============================================  ============================================================