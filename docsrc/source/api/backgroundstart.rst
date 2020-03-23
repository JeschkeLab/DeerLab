.. highlight:: matlab
.. _backgroundstart:

***********************
:mod:`backgroundstart`
***********************

Computes the optimal point for starting the background fit

-----------------------------


Syntax
=========================================

.. code-block:: matlab

    [t0,pos] = backgroundstart(V,t,model)
    [t0,pos] = backgroundstart(V,t,model,'Property',Value)

Parameters
    *   ``V`` - Signal (*N*-element array)
    *   ``t`` - Time axis (*N*-element array)
    *   ``model`` - Background model (function handle)
Returns
    *   ``t0`` - Optimal fitting start time (scalar)
    *   ``pos`` - Optimal fitting start index (scalar)

-----------------------------


Description
=========================================

.. code-block:: matlab

    [t0,pos] = backgroundstart(V,t,model)

Returns the optimal start time ``t0`` and corresponding array index ``pos`` at which to start fitting the background function corresponding to the model given by ``model``. The pre-defined background models defined by the ``model`` argument are the following:

* ``@bg_exp`` - See :ref:`bg_exp` for details.

* ``@bg_strexp`` -  See :ref:`bg_strexp` for details.

* ``@bg_prodstrexp`` -  See :ref:`bg_prodstrexp` for details.

* ``@bg_sumstrexp`` -  See :ref:`bg_sumstrexp` for details.

*  ``@bg_poly1``, ``@bg_poly2``, ``@bg_poly3`` - See :ref:`bg_poly1`, :ref:`bg_poly2` and :ref:`bg_poly3`  for details.

.. image:: ../images/backgroundstart1.svg

-----------------------------


Optional Arguments
=========================================

Optional arguments can be specified by parameter/value pairs. All property names are case insensitive and the property-value pairs can be passed in any order after the required input arguments have been passed..

.. code-block:: matlab

    [t0,pos]  = backgroundstart(args,'Property1',Value1,'Property2',Value2,...)


- ``'SearchStart'`` - Relative Search Start
    Time (in microseconds) at which the background search starts

    *Default:* ``0.1*max(t)``

    *Example:*

		.. code-block:: matlab

			P = backgroundstart(args,'RelSearchStart',1.5)

- ``'SearchEnd'`` - Relative Search End
    Time (in microseconds) at which the background search ends

    *Default:* ``0.6*max(t)``

    *Example:*

		.. code-block:: matlab

			P = backgroundstart(args,'RelSearchEnd',4)

- ``'EndCutOff'`` - Signal cutoff
    Time (in microseconds) after which the signal is no longer used for fitting By default the whole signal is evaluated. 

    *Default:* ``max(t)``

    *Example:*

		.. code-block:: matlab

			P = backgroundstart(args,'EndCutOff',6.5)

- For further property-value pair options see :ref:`fitbackground`.

