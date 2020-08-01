Loading and preprocessing
=========================================

DeerLab provide a wide range of functionality to analyze experimental dipolar EPR data.

Loading
------------------------------------------

To load an experimental time-domain trace, use the function ``deerload``. It supports a variety of spectrometer formats, and it can read in both 1D and 2D data, real- or complex-valued. Here is an example:


.. code-block:: python
   
    t,V,_ = deerload('myfile.DSC')    # load DEER trace and time axis

``deerload`` attempts to return the the time vector ``t`` in correct units (microseconds), but it might not be able to do so for all file formats.


Preprocessing
------------------------------------------
Before fitting, experimental DEER data must be preprocessed. Three steps are usually required:

-  Since the dipolar signal is usually complex, the first step if to perform a phase correction which will minimize the imaginary component / maximize the real component.

- In common commercial spectrometers, the time-axes are measured in absolute values. This step aims to optimally determine and correct for the zero-time of the time axis.

- The intensity of the signal is fiven in some arbitrary units (usually some kind voltage). All functions in DeerAnalysis are agnostic with respect to this scaling and do not require the dipolar signal to be scaled such that V(0)=1. This scale is optimized by all fit functions in DeerLab, yet it can still be corrected by pre-processing.

.. code-block:: python

    V = correctphase(V)           # phase correction
    t = correctzerotime(V,t)      # zero-time adjustment

All the ``correct*`` functions determine the corrections via an optimization approach. If that fails, you can provide an explicit phase, time shift, and scale as additional argument.


Noise levels
------------------------------------------

It can be useful to estimate the level of noise present in the data before actually fitting the data. Use the function ``noiselevel`` for this.

.. code-block:: python

   sigma = noiselevel(V)

This returns an estimate of the noise standard deviation in the signal ``V``. If the signal is 1D, it estimates the noise level from the imaginary part, and if the signal is 2D (containing many scans), it gets from the standard deviation across scans.

