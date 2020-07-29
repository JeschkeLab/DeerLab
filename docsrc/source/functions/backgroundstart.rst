.. highlight:: matlab
.. _backgroundstart:

***********************
:mod:`backgroundstart`
***********************

Computes the optimal point for starting the background fit

.. warning:: This is a legacy function. Its use is not recommended for routine or accurate data analysis.


-----------------------------


Syntax
=========================================

.. code-block:: matlab

    [t0,pos] = backgroundstart(V,t,bgmodel)
    [t0,pos] = backgroundstart(___,'Property',Value)

Parameters
    *   ``V`` - Signal (*N*-element array)
    *   ``t`` - Time axis (*N*-element array)
    *   ``bgmodel`` - Background model (function handle)
Returns
    *   ``t0`` - Optimal fitting start time (scalar)
    *   ``pos`` - Optimal fitting start index (scalar)

-----------------------------


Description
=========================================

.. code-block:: matlab

    [t0,pos] = backgroundstart(V,t,bgmodel)

Returns the optimal start time ``t0`` and corresponding array index ``pos`` at which to start fitting the background function corresponding to the background model given by ``model`` (such as @bg_hom3d etc - see full list of :ref:`background models<modelsref_bg>`). 

-----------------------------


Additional Settings
=========================================

Additional settings can be specified via name-value pairs. All property names are case insensitive and the property-value pairs can be passed in any order after the required input arguments have been passed.

.. code-block:: matlab

    [t0,pos]  = backgroundstart(___,'Property1',Value1,'Property2',Value2,___)


- ``'SearchStart'`` - Relative Search Start
    Time (in microseconds) at which the background search starts

    *Default:* ``0.1*max(t)``

    *Example:*

		.. code-block:: matlab

			[t0,pos] = backgroundstart(___,'RelSearchStart',1.5)

- ``'SearchEnd'`` - Relative Search End
    Time (in microseconds) at which the background search ends

    *Default:* ``0.6*max(t)``

    *Example:*

		.. code-block:: matlab

			[t0,pos] = backgroundstart(___,'RelSearchEnd',4)

- ``'EndCutOff'`` - Signal cutoff
    Time (in microseconds) after which the signal is no longer used for fitting By default the whole signal is evaluated. 

    *Default:* ``max(t)``

    *Example:*

		.. code-block:: matlab

			[t0,pos] = backgroundstart(___,'EndCutOff',6.5)

- For further name-value pair options see :ref:`fitbackground`.

