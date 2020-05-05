.. highlight:: matlab
.. _deerload:

***********************
:mod:`deerload`
***********************

Load experimental DEER data from file


-----------------------------


Syntax
=========================================

.. code-block:: matlab

    V = deerload(filename)
    [t,V] = deerload(filename)
    [t,V,pars,file] = deerload(filename)
    [t,V,pars,file] = deerload


Parameters
    *   ``filename`` -- Name of data file (string)
Returns
    *   ``t`` -- Time axis (*N*-element array)
    *   ``V`` -- Experimental signal (*N*-element array)
    *   ``pars`` -- Parameter file entries (struct)
    *   ``file`` -- Full path to data file(string)




-----------------------------


Description
=========================================

.. code-block:: matlab

    V = deerload(filename)
    [t,V] = deerload(filename)

Read spectral data from a file specified in the string ``filename`` and returns the time-axis vector ``t`` and data vector ``V``.

.. Important::
   Most commercial spectrometers save their data in nanoseconds. Since the required time unit in DeerLab is microseconds, it is important to check the values of ``t`` returned by ``deerload`` and convert to microseconds if necessary.

-----------------------------


.. code-block:: matlab

    [t,V,pars] = deerload(filename)

The structure ``pars`` contains entries from the parameter file, if present.


-----------------------------


.. code-block:: matlab

    [t,V] = deerload


If ``filename`` is a directory, a file browser is displayed. If ``filename`` is omitted, the current directory is used as default. ``deerload`` returns the name of the loaded file (including its path) as fourth parameter ``file``.


