.. highlight:: matlab
.. _noiselevel:

*********************
:mod:`noiselevel`
*********************

Estimate the noise standard deviation on a signal.

-----------------------------


Syntax
=========================================

.. code-block:: matlab

   sigma = noiselevel(V2D)
   sigma = noiselevel(Vc)
   sigma = noiselevel(V)
   sigma = noiselevel(V,filter)
   sigma = noiselevel(V,Vref)

Parameters
    *   ``V2D`` - 2D-dataset of scans of a dipolar signal (*NxM*-element matrix)
    *   ``Vc`` - 1D complex valued dipolar signal (*N*-element array)
    *   ``V`` - 1D dipolar signal (*N*-element array)
    *   ``filter`` - Filtering method (string)
    *   ``Vref`` - Reference signal (*N*-element array)


Returns
    *  ``sigma`` - Estimated noise standard deviation (scalar)

-----------------------------


Description
=========================================

.. code-block:: matlab

   sigma = noiselevel(V2D)

If a 2D-dataset ``V2D`` of different scans for a signal is provided, the noise standard deviation is estimated from the deviations between scans. The second dimension of ``V2D`` must contain the different scans. The function returns the standard deviation of the averaged signal not of the individual scans.


-----------------------------

.. code-block:: matlab

   sigma = noiselevel(Vc)

If a complex-valued signal ``Vc`` is provided, the imaginary component is minimzed via phase correction and the noise standard deviation is estimated from the phase-corrected imaginary component.


-----------------------------

.. code-block:: matlab

   sigma = noiselevel(V)
   sigma = noiselevel(V,filter)

If a real-valued signal ``V`` is provided, the noise standard deviation is estimated from the deviation obtained via application of a moving mean filter. The nature of the filter can specified by setting ``filter`` to one of the following methods: 


	*  ``'movmean'`` - Moving mean filter (default)
	*  ``'savgol'`` - Savitzky-Golay filter
	

-----------------------------

.. code-block:: matlab

   sigma = noiselevel(V,Vref)

If a reference model signal ``Vref`` is given, the noise level is estimated from the difference between both signals.