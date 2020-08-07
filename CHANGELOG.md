
-------------------------------

Version 0.10.0 - August 2020
-----------------------------

As of this version, DeerLab is written in Python in contrast to older versions based on MATLAB.

Deprecated functions
    The following functions have been deprecated due to limited usability or due to functionality overlap with other DeerLab functions: ``aptkernel``, ``backgroundstart``, ``fitbackground``, ``paramodel``, and ``time2freq``. 

Overall changes
    * All fit functions now return a single ``FitResult`` output which will contain all results. 

    * All functions are now compatible with non-uniformly increasing distance axes. 

    * All fit functions are completely agnostic with respect of the abolute values of the signal amplitude. This is automatically fitted by all function and return as part of the results.

    * Uncertainty quantification for all fit functions is returned as a ``UncertQuant`` object from which confidence intervals, parameter distributions, etc. can be generated generalizing the uncertainty interface for all DeerLab. Uncertainty can now be propagated to arbitrary functions.

Specific changes
    * ``fitparamodel``: the functionality has been streamlined. Function now fits arbitrary parametric models using non-linear leas-squares without consideration of whether it is a time-domain or distance-domain model. The models no longer need to take two inputs (axis+parameters) and now only tk the parameters as input. 

    * ``fitregmodel``: goodness-of-fit estimators are now computed using the proper estimation the degrees of freedom.

    * ``fitmultimodel``: added internal measures to avoid situations where one or several components are suppressed by fitting zero-amplitudes making the method more stable. 

    * ``uqst``: the uncertainty distributions of the parameters are now returned as properly normalized probability density functions.

    * ``fftspec``: frequency axis construction has been corrected.

    * ``regoperator``: now calculates the numerically exact finite-difference matrix using Fornberg's method.

    * ``correctphase``: now can handle 2D-datasets.

-------------------------------