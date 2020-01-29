.. highlight:: matlab
.. _deerload:

***********************
:mod:`deerload`
***********************

Load experimental spectrometer data


-----------------------------


Syntax
=========================================

.. code-block:: matlab

    V = deerload(filename)
    [t,V] = deerload(filename)
    [t,V,pars] = deerload(filename)
    [t,V,pars,file] = deerload(filename)
    [t,V] = deerload


Parameters
    *   ``filename`` - Name of data file (string)
    *   ``t`` - Time axis (*N*-element array)
Returns
    *   ``t`` - Time axis (*N*-element array)
    *   ``V`` - Experimental signal (*N*-element array)
    *   ``pars`` - Parameter file entries (struct)
    *   ``file`` - Full path to data file(string)




-----------------------------


Description
=========================================

.. code-block:: matlab

    V = deerload(filename)
    [t,V] = deerload(filename)

Read spectral data from a file specified in the string ``filename`` into the arrays ``t`` (time-axis in microseconds) and ``V`` (data).


-----------------------------


.. code-block:: matlab

    [t,V,pars] = deerload(filename)

The structure ``pars`` contains entries from the parameter file, if present.


-----------------------------


.. code-block:: matlab

    [t,V] = deerload


If ``filename`` is a directory, a file browser is displayed. If ``filename`` is omitted, the current directory is used as default. ``deerload`` returns the name of the loaded file (including its path) as fourth parameter ``file``.


