
Release v0.13.0 - April 2021
---------------------------------

#### New features 

- DeerLab now supports Python 3.9.
- The function ``fitsignal`` has been re-named to ``fitmodel`` for correctness and consistency with other functions ([#102](https://github.com/JeschkeLab/DeerLab/pull/102)).
- Added new experiment models for RIDME on systems with one to seven harmonic pathways (S=1/2 to S=7/2) to include all higher harmonics (overtones) ([#79](https://github.com/JeschkeLab/DeerLab/pull/79)). 
- Bootstrapping is now embedded into ``fitmodel`` to automatically bootstrap all output quantities without the need to write additional script lines ([#55](https://github.com/JeschkeLab/DeerLab/issues/55)). In ``fitmodel`` a new option ``uq`` allows to switch between covariance or bootstrapping uncertainty quantification ([#88](https://github.com/JeschkeLab/DeerLab/pull/88)). 
- The function ``fitmodel`` now returns ``Vmod`` and ``Vunmod``, the modulated and unmodulated contributionsto the fitted dipolar signal, respectively, along their uncertainties as additional outputs ([#78](https://github.com/JeschkeLab/DeerLab/pull/78)).
- Implemented several initialization strategies in ``fitmultimodel`` for multi-model components ([#67](https://github.com/JeschkeLab/DeerLab/pull/67)). Three different new strategies ``'spread'``, ``'split'`` and ``'merge'`` will initialize the parameter values of the N-component fit based on the results of the N-1/N+1 component fit to improve quality of results and speed.  
- Added contribution guidelines to the documentation and automated list of DeerLab contributors. 
- The function ``snlls`` now accepts additional custom penalties to include in the optimization ([#76](https://github.com/JeschkeLab/DeerLab/issues/76), [#108](https://github.com/JeschkeLab/DeerLab/pull/112)).
- All fit functions now return the fit of the data along its uncertainty automatically as part of the ``FitResult`` object ([#130](https://github.com/JeschkeLab/DeerLab/issues/130), [#134](https://github.com/JeschkeLab/DeerLab/pull/134)).

#### Overall changes

- The performance of all fit functions has been considerably accelerated by removing call overheads in built-in DeerLab models ([#100](https://github.com/JeschkeLab/DeerLab/issues/100), [#101](https://github.com/JeschkeLab/DeerLab/pull/101), [#143](https://github.com/JeschkeLab/DeerLab/pull/143)).
- Improved robustness of the installer ([#65](https://github.com/JeschkeLab/DeerLab/pull/65)):
    - The installer no longer assumes the alias ``pip`` to be setup on the system. 
    - The installation will now handle cases when system-wide privileges are not available ([#52](https://github.com/JeschkeLab/DeerLab/issues/52)).
    - Improved robustness of the installation in Windows systems to avoid missing DLL errors ([#64](https://github.com/JeschkeLab/DeerLab/issues/64)).
    - The installer will now get the latest Numpy/Scipy releases in Windows systems available at the [Gohlke repository](https://www.lfd.uci.edu/~gohlke/pythonlibs/). 
- Adapted piece of code leading to a ``VisibleDeprecationWarning`` visible during execution of certain DeerLab functions.
- Improved interface of built-in plots in ``FitResult.plot()``. The method now returns a Matplotlib figure object (`matplotlib.figure.Figure`) instead of an axes object (``matplotlib.axes._subplots.AxesSubplot``) which can be modified more freely to adjust graphical elements ([#85](https://github.com/JeschkeLab/DeerLab/issues/85)). The method now takes an optional keyword ``FitResult.plot(show=True\False)`` to enable/disable rendering of the graphics upon calling the method ([#87](https://github.com/JeschkeLab/DeerLab/pull/87)).
- The fit objective values returned in ``FitResult.cost`` are now correct (previous versions had an erroneous 1/2 factor) ([#80](https://github.com/JeschkeLab/DeerLab/issues/80)). The value is now returned as a scalar value instead of a single-element list ([#81](https://github.com/JeschkeLab/DeerLab/issues/81)).
- Removed the re-normalization conventions ``K(t=0,r)=1`` and ``B(t=0)=1`` and associated options ``renormalize`` and ``renormpaths`` in the ``dipolarkernel`` and ``dipolarbackground`` functions ([#99](https://github.com/JeschkeLab/DeerLab/pull/99)) to avoid identifiability issues between dipolar pathway amplitudes and signal scales during fitting ([#76](https://github.com/JeschkeLab/DeerLab/issues/76)). 
- The fit convergence criteria ``tol`` (objective function tolerance) and ``maxiter`` (iteration limit) are now exposed as keyword argument in all fit functions ([#111](https://github.com/JeschkeLab/DeerLab/issues/111), [#112](https://github.com/JeschkeLab/DeerLab/pull/112)). 
- Multiple improvements and corrections to the documentation ([#95](https://github.com/JeschkeLab/DeerLab/pull/95), [#96](https://github.com/JeschkeLab/DeerLab/pull/96), [#104](https://github.com/JeschkeLab/DeerLab/pull/104), [#106](https://github.com/JeschkeLab/DeerLab/pull/106), [#107](https://github.com/JeschkeLab/DeerLab/pull/107), [#115](https://github.com/JeschkeLab/DeerLab/pull/115), [#122](https://github.com/JeschkeLab/DeerLab/pull/122))
- Corrections in the metadata of multiple ``dd_models``. The key ``Parameters`` of some models contained the wrong names.
- The metadata of the built-in models is now accessible and manipulable via function attributes (e.g. ``dd_gauss.parameters``) rather than trought a returned dictionary (e.g. ``dd_gauss()['Parameters']``) ([#143](https://github.com/JeschkeLab/DeerLab/pull/143)).
- The keyword argument to request uncertainty quantification has been unified across all fitting functions. It is now ``uq`` ([#120](https://github.com/JeschkeLab/DeerLab/pull/120)).
- The ``UncertQuant`` class has been renamed into ``UQResult`` ([#123](https://github.com/JeschkeLab/DeerLab/pull/123)).
- Uncertainty quantification is now tested numerically against an external package (``lmfit``) to ensue quality and accuracy ([#121](https://github.com/JeschkeLab/DeerLab/pull/121)).

#### Specific changes
- ``deerload``: 
    - Fixed behaviour of the function when loading certain 2D-datasets in the BES3T format ([#82](https://github.com/JeschkeLab/DeerLab/issues/82), [#83](https://github.com/JeschkeLab/DeerLab/pull/83)).
    - In 2D-datasets, the abscissas are now returned as a list of abscissas instead of a single 2D-matrix ([#83](https://github.com/JeschkeLab/DeerLab/pull/83))). 
- ``fitmodel``:
    - Corrected the scaling behaviour of all outputs related to components of the dipolar signal to match the scaling of the original experimental data ([#78](https://github.com/JeschkeLab/DeerLab/pull/78)). 
    - The built-in plot method ``FitResult.plot()`` now plots the unmodulated component fit as well with its uncertainty ([#78](https://github.com/JeschkeLab/DeerLab/pull/78)).
    - Extended information included in the verbose summary ([#78](https://github.com/JeschkeLab/DeerLab/pull/78)). 
    - Simplified the interface for defining initial values and boundaries of parameters in ``fitsignal`` ([#71](https://github.com/JeschkeLab/DeerLab/pull/71)). Now instead of defining, e.g., `fitsignal(..., lb = [[],[50],[0.2, 0.5]])` one can define the individual vales/boundaries ``fitsignal(..., bg_lb = 50, ex_lb = [0.2,0.5])`` by using the new keywords. 
    - Removed the keyword argument ``uqanalysis=True/False``. The uncertainty quantification can now be disabled via the new keyword ``uq=None`` ([#98](https://github.com/JeschkeLab/DeerLab/pull/98)).
    - Corrected the behaviour of built-in start values when manually specifying boundaries ([#73](https://github.com/JeschkeLab/DeerLab/pull/73)). If the built-in start values are outside of the user-specified boundaries the program will now automatically set the start values in the middle of the boundaries to avoid errors ([#72](https://github.com/JeschkeLab/DeerLab/issues/72)).
    - Implemented the constraint ``Lam0+sum(lam)<=1`` to ensure the structural-identifiability of ``Lam0`` and ``V0`` during SNLLS optimization of experiment models with more than one modulated dipolar pathway (i.e. does not affect ``ex_4pdeer``) ([#76](https://github.com/JeschkeLab/DeerLab/issues/76),[#108](https://github.com/JeschkeLab/DeerLab/pull/108)).
    - Removed the optional keyword argument ``regtype`` ([#137](https://github.com/JeschkeLab/DeerLab/pull/137)).
- ``fitregmodel``:
    - Corrected the behaviour of the uncertainty quantification when disabling the non-negativity constraint ([#121](https://github.com/JeschkeLab/DeerLab/pull/121)).
- ``fitparamodel``: 
    - Made ``par0`` a positional argument instead of an optional keyword ([#70](https://github.com/JeschkeLab/DeerLab/issues/70)). to avoid errors when not defined ([#69](https://github.com/JeschkeLab/DeerLab/issues/69)).
    - Keyword argument ``rescale`` has been renamed to ``fitscale`` ([#128](https://github.com/JeschkeLab/DeerLab/issues/128),[#129](https://github.com/JeschkeLab/DeerLab/pull/129)).
- ``snlls``:
    - Corrected bug that was leading to the smoothness penalty being accounted for twice in the least-squares residual during optimization ([#103](https://github.com/JeschkeLab/DeerLab/issues/103)).
    - Now returns the uncertainty quantification of linear and nonlinear parts as separate objects ``nonlinUncert`` and ``linUncert`` ([#108](https://github.com/JeschkeLab/DeerLab/pull/108)).
    - Improved the covariance-based uncertainty analysis by including correlations between linear and non-linear parameters ([#108](https://github.com/JeschkeLab/DeerLab/pull/108)).
    - Improved the behavior of signal scale determination ([#108](https://github.com/JeschkeLab/DeerLab/pull/108)).
    - Enabled prescaling of the data to avoid scaling issues during uncertainty quantification ([#132](https://github.com/JeschkeLab/DeerLab/issue/132), [#133](https://github.com/JeschkeLab/DeerLab/pull/133)).
    - Corrected the behaviour of the uncertainty quantification when disabling the regularization penalty ([#121](https://github.com/JeschkeLab/DeerLab/pull/121)).
- ``diststats``: 
    - Now compatible with non-uniformly defined distance distributions ([#92](https://github.com/JeschkeLab/DeerLab/issues/92), [#94](https://github.com/JeschkeLab/DeerLab/pull/94)). 
    - Added internal validation step to avoid non-sensical results when confounding the syntax ([#91](https://github.com/JeschkeLab/DeerLab/pull/91)).
- ``dipolarkernel``: 
    - Now allows defining pathways without unmodulated components.
    - All optional keyword arguments can only be passed as named and not positional arguments ([#138](https://github.com/JeschkeLab/DeerLab/pull/138)). 
    - The keyword ``pathways`` now only takes lists of pathways and not modulation depth parameters. A new separate keyword ``mod`` takes the modulation depth parameter for the simplified 4-pulse DEER kernel ([#118](https://github.com/JeschkeLab/DeerLab/issues/118), [#138](https://github.com/JeschkeLab/DeerLab/pull/138)).
    - Renmaed the background argument keyword ``B`` into ``bg`` ([#138](https://github.com/JeschkeLab/DeerLab/pull/138)).
- ``regparamrange``:
    - Implemented new CSD algorithm to avoid LAPACK library crashes encountered when using multiple DeerLab functions calling ``regparamrange`` internally ([#68](https://github.com/JeschkeLab/DeerLab/pull/68)).
- ``correctphase``: 
    - Implement new keyword ``phase`` to select the criterion for optimizing the phase for correction ([#114](https://github.com/JeschkeLab/DeerLab/issues/114), [#131](https://github.com/JeschkeLab/DeerLab/pull0/131)).
    - Deprecated imaginary offset fitting ([#131](https://github.com/JeschkeLab/DeerLab/pull0/131)). 
    - Deprecated manual phase correction. Manual correction can be done by the user and is now described in the beginner's guide ([#131](https://github.com/JeschkeLab/DeerLab/pull0/131)). 
    
-------------------------------

Release v0.12.2 - October 2020
---------------------------------

#### Hotfix

* ``regparamrange``: The exception handling introduced in the previous release was still too specific. The function kept crashing due to SVD non-convergence errors during the GSVD. This has been fixed and the error will not lead to a crash. ([#42](https://github.com/JeschkeLab/DeerLab/issues/42)).   

#### Overall changes

* Fit functions using the ``multistart`` option are now fully deterministic. The functions was using now a random generator to define the different start points, this is now deterministic. 

* Documentation UI has been re-designed for a more confortable reading. Minor errors and outdated information have been corrected throughout. Expanded reference documentation of several functions for better understanding. 

#### Specific changes

* ``dd_skewgauss``: corrected an error in the implementation that was leading to wrong distributions ([#61](https://github.com/JeschkeLab/DeerLab/issues/61)).  

* ``dd_models``, ``ex_models``: Adapted numerical boundaries and start values of some built-in models to reflect better the physical reality. Afected models: ``dd_skewgauss``, ``dd_triangle``, ``dd_gengauss``, ``ex_5pdeer``, ``ex_ovl4pdeer``. 

-------------------------------

Release v0.12.1 - October 2020
---------------------------------

#### Overall changes

* The calculation of the Jacobian for covariance-based uncertainty analysis has been simplified providing a significant boost in performance for all fit functions ([#55](https://github.com/JeschkeLab/DeerLab/pull/55)). 

* The Jacobian computation is more robust, now taking into consideration parameter boundaries ([#58](https://github.com/JeschkeLab/DeerLab/pull/58)). This solves errors such as the ones reported in [#50](https://github.com/JeschkeLab/DeerLab/issues/50).

* Broken examples in the documentation have been corrected ([#57](https://github.com/JeschkeLab/DeerLab/pull/57)).

* When requesting attributes or method of a UncertQuant object under disabled uncertainty analysis (``uqanalysis=False``) now it will prompt an explanatory error instead of just crashing ([#56](https://github.com/JeschkeLab/DeerLab/issues/56)). 

#### Specific changes

* ``fitsignal``: corected the behaviour of the scaling output (``fit.scale``). Now all fitted dipolar signals (``fit.V``) have the same scaling as the input signal ([#53](https://github.com/JeschkeLab/DeerLab/issues/53)). 

* ``regparamrange``: relaxed the exception handling to catch errors occuring under certain conditions. The function seems to crash due to LAPACK or SVD non-convergence errors during the GSVD, now these are catched and the alpha-range is estimated using simple SVD as an approximation. This function might be deprecated in a future release ([#42](https://github.com/JeschkeLab/DeerLab/issues/42)).   

-------------------------------


Release v0.12.0 - October 2020
---------------------------------

#### Overall changes

* The regularization operator matrices ``regoperator`` now include the edges of the distribution ([#38](https://github.com/JeschkeLab/DeerLab/pull/38)). Now the smoothness penalty is imposed on the distribution edges avoiding the accumulation of distribution mass at the edges of ``r``. 

* The interface for defining dipolar pathways has been simplified ([#41](https://github.com/JeschkeLab/DeerLab/pull/41)). For example, a signal with two dipolar pathways had to be defined as ``pathways = [[Lam0,np.nan], [lam1,T0]]``. Now the unmodulated pathway must be defined by its amplitude and does not accept the use of ``np.nan``, e.g. ``pathways = [Lam0, [lam1,T0]]``.

* The project version control has been switched from the Git-flow to the GitHub-flow design. The default branch has been switched from ``master`` to ``main``, which is now always production-ready. All new contributions are merged into ``main`` exclusively by pull requests.

* The dependency on the ``lambda`` parameter has been removed from all phenomenological background models, and kept only for physical models ([#43](https://github.com/JeschkeLab/DeerLab/pull/43)). Their interface with ``dipolarbackground`` and ``dipolarkernel`` have been updated accordingly. 

#### New features: 

* ``diststats``: New function to calculate different statistical quantities of the distance distribution and their corresponding uncertainties ([#37](https://github.com/JeschkeLab/DeerLab/pull/37)).

* ``bootan``: Introduced the option ``cores`` to parallelize the bootstrapping using multiple CPUs ([#35](https://github.com/JeschkeLab/DeerLab/pull/35)). 

#### Specific changes

* ``bg_homfractal``: 
    -  Corrected behavior of the model. For ``d=3`` the model returned wrong values, and for ``d~=3`` the model resulted in an error.

* ``UncertQuant``: 
    - Fixed bug when propagating uncertainty to scalar functions.

* ``deerload``: 
    - Fixed UTF-8 error when loading certain spectrometer files in MacOS ([#30](https://github.com/JeschkeLab/DeerLab/pull/30))

* ``fitsignal``:
    - The fitted scale of the signal is now properly calculated when fitting fully parametric signals. 
    - Fixed error occuring when fitting a dipolar evolution function with a non-parametric distribution.

* ``selregparam``:
    - Fixed bug occuring when requesting the ``lc`` or ``lr`` selection methods.

* ``regparamange``: 
    - An error occuring at the BLAS/LAPLACK error ocurring under certain conditions in MacOS and Ubuntu is now handled to avoid a crash. 

-------------------------------


Release v0.11.0 - September 2020
---------------------------------

#### Overall changes

* All Gauss models (``dd_gauss``,etc.) now use the standard deviation ``sigma`` instead of the FWHM as the width parameter for consistency with other method such as Rice distributions ([#19](https://github.com/JeschkeLab/DeerLab/pull/19)).

* All hard-wired random seeds have been removed. 

* A new method ``plot()`` has been added to the ``FitResult`` class returned by all fit functions. This will create a basic plot of the fit results ([#7](https://github.com/JeschkeLab/DeerLab/pull/7)).

#### Specific changes
* ``snlls``: 
    - Renamed option ``penalty`` as ``reg`` and improved its interface ([#13](https://github.com/JeschkeLab/DeerLab/pull/13)).
    - The regularization parameter of the optimal solution is returned now ([#20](https://github.com/JeschkeLab/DeerLab/pull/20)).

* ``whitegaussnoise``: Added a ``seed`` option to select static noise realizations.

* ``correctzerotime``: 
    - Fixed bug when zero-time is at start/end of array ([#24](https://github.com/JeschkeLab/DeerLab/pull/24)).
    - Function no longer rescales the experimental data passed on to the function. 

* ``fitsignal``:  
    - The regularization parameter of the optimal solution is returned now ([#20](https://github.com/JeschkeLab/DeerLab/pull/20)).
    - Bug fixed when fitting dipolar evolution functions (no background and no experiment models) with a parametric distance distribution. 

* ``fitmultimodel``: Start points are now spread over constrained parameter space grid instead of being randomble initiated ([#22](https://github.com/JeschkeLab/DeerLab/pull/22)).

* ``deerload``: 
    - Now returns the time axis in microseconds instead of nanoseconds ([#21](https://github.com/JeschkeLab/DeerLab/pull/21)).
    - The bug appearing when loading certain BES3T files has been fixed ([#14](https://github.com/JeschkeLab/DeerLab/pull/14)).

* ``fitregmodel``: Now returns the fitted dipolar signal in the ``FitResult`` output

* ``correctscale``: The parameter fit ranges have been adjusted.


-------------------------------

Release v0.10.0 - August 2020
-----------------------------

As of this version, DeerLab is based on Python in contrast to older versions based on MATLAB found [here](https://github.com/JeschkeLab/DeerLab-Matlab).

#### Deprecated functions
The following functions have been deprecated due to limited usability or due to functionality overlap with other DeerLab functions: ``aptkernel``, ``backgroundstart``, ``fitbackground``, ``paramodel``, and ``time2freq``. 

#### Overall changes
* All fit functions now return a single ``FitResult`` output which will contain all results. 

* All functions are now compatible with non-uniformly increasing distance axes. 

* All fit functions are completely agnostic with respect of the abolute values of the signal amplitude. This is automatically fitted by all function and return as part of the results.

* Uncertainty quantification for all fit functions is returned as a ``UncertQuant`` object from which confidence intervals, parameter distributions, etc. can be generated generalizing the uncertainty interface for all DeerLab. Uncertainty can now be propagated to arbitrary functions.

#### Specific changes
* ``fitparamodel``: the functionality has been streamlined. Function now fits arbitrary parametric models using non-linear leas-squares without consideration of whether it is a time-domain or distance-domain model. The models no longer need to take two inputs (axis+parameters) and now only tk the parameters as input. 

* ``fitregmodel``: goodness-of-fit estimators are now computed using the proper estimation the degrees of freedom.

* ``fitmultimodel``: added internal measures to avoid situations where one or several components are suppressed by fitting zero-amplitudes making the method more stable. 

* ``uqst``: the uncertainty distributions of the parameters are now returned as properly normalized probability density functions.

* ``fftspec``: frequency axis construction has been corrected.

* ``regoperator``: now calculates the numerically exact finite-difference matrix using Fornberg's method.

* ``correctphase``: now can handle 2D-datasets.

