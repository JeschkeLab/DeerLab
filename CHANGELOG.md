
-------------------------------

Patch v0.12.1 - October 2020
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

* The interface for defining dipolar pathways has been simplified ([#41](https://github.com/JeschkeLab/DeerLab/pull/41)). For example, a signal with two dipolar pathways had to be defined as ``pathways = [[Lam0,np.nan], [lam1,T0]]``. Now the unmodulated pathway can just be defined by its amplitude and does not require the use of ``np.nan``, e.g. ``pathways = [Lam0, [lam1,T0]]``.

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

-------------------------------