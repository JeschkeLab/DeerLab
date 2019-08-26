.. highlight:: matlab
.. _fitbackground:


***********************
:mod:`fitbackground`
***********************

Fit the background function in a signal

"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
:mod:`[B,param] = fitbackground(S,t,tfit,'model',p)`
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Parameters
    *   **S** - Data to fit (M-array)
    *   **t** - Time axis (N-array)
    *   **tfit** - Time axis to fit (M-array)
    *   **model** - Background model (string)
    *   **p** - Additional model parameters (array)

Returns
    *   **B** - Spectrum (M-array)
    *   **param** - Fitted parameter values (array)
Usage
=========================================

.. code-block:: matlab

    [B,param] = fitbackground(S,t,tfit,'model')

Fits the the paramters ``param`` of the N-point background function ``B``. This is done by fitting the M-point data ``S`` on a M-point axis ``tfit`` using a model given by the string ``'model'``. The background is then extrapolated to the N-point axis ``t``. The pre-defined models in fitbackground defined by the 'model' string argument are the following:

* ``'exponential'`` - exponential function where the decay rate if fitted

* ``'polyexp'`` -  exponetial function by fitting the a linear function to the log of the signal

* ``'fractal'`` - stretched exponential by fitting the decay rate and fractal dimension on the log of the signal

* ``'polynomial'`` - polynomial function of order as given as an input

.. code-block:: matlab

    [B,param] = fitbackground(S,t,tfit,'polynomial',p)

For polynomial function fitting, the polynomial order can be specified as an additional input argument ``p``.

.. code-block:: matlab

    [B,param] = fitbackground(S,t,tfit,@model)

User-defined models can be fittid by passing a function handle instead of a model name. To pass user-defined models, the @model argument must be a function handle to a function accepting two input arguments as follows ``function myModel(t,param), ..., end`` where ``param`` is an array containing the parameter of the model. Example models include the ``@strexp``, ``@sumstrexp``, ``@prodstrexp`` models distibuted in DeerAnalysis2.
