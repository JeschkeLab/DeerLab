Loading and preprocessing
=========================================

DeerLab provide a wide range of functionality to analyze experimental dipolar EPR data.

Loading
------------------------------------------

To load an experimental time-domain trace, use the function ``deerload``. It supports a variety of spectrometer formats, and it can read in both 1D and 2D data, real- or complex-valued. Here is an example:


.. code-block:: python
   
    t,V = dl.deerload('myfile.DSC')    # load DEER trace and time axis

``deerload`` attempts to return the the time vector ``t`` in correct units (microseconds), but it might not be able to do so for all file formats.


Preprocessing
------------------------------------------
Before fitting, experimental DEER data must be preprocessed. Three steps are usually required:

- Phase correction. Since the dipolar signal is usually complex, the first step if to perform a phase correction which minimizes the imaginary component and maximizes the real component.

- Zero-time correction. Experimental time-axes values might be shifted relative to the theoretical ones. This step determines the zero-time of the time axis based on the data and shifts the time axis.

.. code-block:: python

    V = dl.correctphase(V)           # phase correction
    t = dl.correctzerotime(V,t)      # zero-time adjustment

Both ``correct*`` functions determine the corrections via an optimization approach. If that fails, you can provide an explicit phase, time shift, and scale as additional argument.

The signal amplitude is given in arbitrary units. DeerLab functions are agnostic with respect to this scaling and do not require the dipolar signal to be scaled. The scale is automatically optimized by all fit functions in DeerLab.


Noise levels
------------------------------------------

It can be useful to estimate the noise level in the data before fitting the data. Use the function ``noiselevel`` for this.

.. code-block:: python

   sigma = dl.noiselevel(V)

This returns an estimate of the noise standard deviation in the signal ``V``. If the signal is 1D, it estimates the noise level from the imaginary part or by filtering, and if the signal is 2D (containing many scans), it gets it from the standard deviation across scans.
