.. highlight:: matlab
.. _backgroundstart:

***********************
:mod:`backgroundstart`
***********************

Computes the optimal point for starting the background fit

"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
:mod:`[t0,pos] = backgroundstart(V,t,model,...)`
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Parameters
    *   **V** - Signal (N-array)
    *   **t** - Time axis (N-array)
    *   **model** - Background model (function handle)
Returns
    *   **t0** - Optimal fitting start time (scalar)
    *   **pos** - Optimal fitting start index (scalar)

Usage
=========================================

.. code-block:: matlab

    [t0,pos] = backgroundstart(V,t,model)

Returns the optimal start time ``t0`` and corresponding array index ``pos`` at which to start fitting the background function corresponding to the model given by ``model``. The pre-defined background models defined by the ``model`` argument are the following:

* ``@td_exp`` - Exponential function where the decay rate if fitted

* ``@td_strexp`` -  Exponential function by fitting the a linear function to the log of the signal

* ``'fractal'`` - Stretched exponential function by fitting the decay rate and fractal dimension on the log of the signal

*  ``@td_poly1``, ``@td_poly2``, ``@td_poly3`` - Polynomial functions of various orders

.. image:: ./images/backgroundstart1.svg

Optional Arguments
=========================================

Optional arguments can be specified by parameter/value pairs. All property names are case insensitive and the property-value pairs can be passed in any order after the required input arguments have been passed..

.. code-block:: matlab

    P = backgroundstart(args,'Property1',Value1,'Property2',Value2,...)

.. centered:: **Property Names & Descriptions**

RelSearchStart
    Relative position at which the background start search starts.

    *Default:* ``0.1``

    *Example:*

    .. code-block:: matlab

       P = selregparam(args,'RelSearchStart',0.25)

RelSearchEnd
    Relative position at which the background start search stops.

    *Default:* ``0.6``

    *Example:*

    .. code-block:: matlab

       P = selregparam(args,'RelSearchEnd',0.7)

ModelParam
    Parameter values for the background model (if required).

    *Default:* [*empty*]

    *Example:*

    .. code-block:: matlab

       P = selregparam(args,'polynomial','ModelParam',2) %Polynomial 2nd order

For further property-value pair options see :ref:`fitbackground`.

