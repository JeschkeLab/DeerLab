.. _functions:


Functions
---------------------

The following list contains the names of the different function and a brief description of their functionality. The parametric model functions are listed in separate sections.

Modeling
=========================================

This class of functions allows simulation of dipolar signals and their modeling for fitting experimental data.

.. toctree::
    :hidden:
    :glob:
    :maxdepth: 1

    ./functions/dipolarkernel
    ./functions/dipolarbackground
    ./functions/whitenoise
    ./functions/mixmodels

.. rst-class:: func-list

=============================================  ============================================================
Function                                         Description
=============================================  ============================================================
:ref:`dipolarkernel`                           Dipolar kernel constructor
:ref:`dipolarbackground`                       Multi-pathway background constructor
:ref:`whitegaussnoise`                         Gaussian white noise generator
:ref:`mixmodels`                               Parametric model mixer
=============================================  ============================================================


Analysis
=========================================

This class of functions can be used and/or combined to create fitting routines of dipolar data.

.. toctree::
    :hidden:
    :glob:
    :maxdepth: 1

    ./functions/snlls
    ./functions/fitmodel
    ./functions/fitmultimodel
    ./functions/fitparamodel
    ./functions/fitregmodel
    ./functions/regoperator
    ./functions/bootan


.. rst-class:: func-list

=============================================  ============================================================
Function                                         Description
=============================================  ============================================================
:ref:`snlls`                                     Separable non-linear least squares solver
:ref:`fitmodel`                                 Full model signal fitting engine
:ref:`bootan`                                    Bootstrap uncertainty analysis
:ref:`fitmultimodel`                             Multi-component model fitting engine
:ref:`fitparamodel`                              Parametric model fitting engine
:ref:`fitregmodel`                               Regularization fitting engine
:ref:`regoperator`                               Regularization operator constructor
=============================================  ============================================================


Pre-Processing
=========================================

This class of functions provide tools for preparing experimental data for analysis. 

.. toctree::
    :hidden:
    :glob:
    :maxdepth: 1

    ./functions/correctphase
    ./functions/correctzerotime
    ./functions/correctscale

.. rst-class:: func-list

=============================================  ============================================================
Function                                         Description
=============================================  ============================================================
:ref:`correctphase`                              IQ Phase correction
:ref:`correctzerotime`                           Dipolar zero-time correction
=============================================  ============================================================



Model Selection
=========================================

This class of functions helps to find an optimal choice of model or model parameters.

.. toctree::
    :hidden:
    :glob:
    :maxdepth: 1

    ./functions/selregparam
    ./functions/regparamrange


.. rst-class:: func-list

=============================================  ============================================================
Function                                         Description
=============================================  ============================================================
:ref:`selregparam`                              Regularization parameter selector
:ref:`regparamrange`                            Regularization parameter range selector
=============================================  ============================================================


Utilities
=========================================

This class of functions provides several tools for quick commands typically required in data processing.


.. toctree::
    :hidden:
    :glob:
    :maxdepth: 1

    ./functions/deerload
    ./functions/noiselevel
    ./functions/fftspec
    ./functions/diststats



.. rst-class:: func-list

=============================================  ============================================================
Function                                         Description
=============================================  ============================================================
:ref:`deerload`                                  Spectrometer data loader
:ref:`noiselevel`                                Noise level estimator
:ref:`fftspec`                                   Fast-Fourier transform spectrum
:ref:`diststats`                                 Distance distribution statistics
=============================================  ============================================================

---------------------------

Legacy Functions
=========================================

.. toctree::
    :hidden:
    :glob:
    :maxdepth: 1

    ./functions/time2dist
    ./functions/correctscale

This group of functions provides tools for reproducing analysis methods or workflows encountered in older software. These functions are not recommended for routine data analysis.

.. rst-class:: func-list

=============================================  ============================================================
Function                                         Description
=============================================  ============================================================
:ref:`time2dist`                                 Heuristic time-to-distance axis converter
:ref:`correctscale`                              Dipolar signal amplitude rescaling
=============================================  ============================================================
