.. highlight:: matlab
.. _dipolarsignal:

*********************
:mod:`dipolarsignal`
*********************

Generation of dipolar signals from distance distributions

-----------------------------



Syntax
=========================================

.. code-block:: matlab

    V = dipolarsignal(t,r,P)
    V = dipolarsignal(t,r)
    V = dipolarsignal(t,r,P,'Property',Value)


Parameters
    *   ``t`` - Time axis (*N*-element array)
    *   ``r`` - Distance axis (*M*-element array)
    *   ``P`` - Distance distribution (*M*-element array)

Returns
    *   ``V`` - Dipolar signal (N-array)

-----------------------------



Description
=========================================

.. code-block:: matlab

    V = dipolarsignal(t,r,P)

Generates the noiseless form factor ``V`` on a time axis ``t`` from the distance distribution ``P`` in a distance axis ``r``.


-----------------------------


.. code-block:: matlab

    V = dipolarsignal(t,r)

If no distribution ``P`` is provided, then the dipolar signal corresponding to a single distance ``r`` is computed.

-----------------------------



Optional Arguments
=========================================
Optional arguments can be specified by parameter/value pairs. All property names are case insensitive and the property-value pairs can be passed in any order after the required input arguments have been passed..

.. code-block:: matlab

    V = dipolarsignal(args,'Property1',Value1,'Property2',Value2,...)


- ``'ModDepth'`` Modulation depth
    Modulation depth of the form factor, must be passed as a scalar value between 0 and 1.

    *Default:* ``1``

    *Example:*

		.. code-block:: matlab

			V = dipolarsignal(args,'ModDepth',0.4)


- ``'Background'`` - Background function
    Background function passed as a *N*-element array.

    *Default:* [*empty*]]

    *Example:*

		.. code-block:: matlab

			V = dipolarsignal(args,'Background',bg_srtexp(t,param))

- ``'NoiseLevel'`` - Level of noise on the signal
    Scalar value containing the desired standard deviation of a Gaussian noise vector 

    *Default:* ``0``

    *Example:*

		.. code-block:: matlab

			V = dipolarsignal(args,'NoiseLevel',0.05)

		.. Important::
			Each call of ``dipolarsignal`` will return a different noise realization. To set the output to a fixed noise realization, the random number generator must be fixed. In MATLAB this can be accomplished by calling ``rng(k)`` where ``k`` is some integer number.


- ``'Overtones'`` - RIDME overtone coefficients
    Array of RIDME overtone coefficients. The coefficients must be normalized, i.e. they must sum to unity.

    *Default:* [*empty*]

    *Example:*

		.. code-block:: matlab

			V = dipolarsignal(args,'Overtones',[0.2 0.5 0.3])

- ``'g'`` - Electron g-value
    Specifies the g-value of the electron spin center used to compute the dipolar frequencies from the given distance axis.

    *Default:* ``2.004602204236924``

    *Example:*

		.. code-block:: matlab

			K = dipolarkernel(args,'g',2.00) %Use experimental g-value

- ``'Scale'`` - Amplitude scale
    Vertical scale to multiply to the output signal

    *Default:* ``1``

    *Example:*

		.. code-block:: matlab

			V = dipolarsignal(args,'Scale', 1e3)

- ``'Phase'`` - IQ phase of the signal
    Scalar-valued phase of the complex-valued signal (in radians).

    *Default:* ``0``

    *Example:*

		.. code-block:: matlab

			V = dipolarsignal(args,'Phase', pi/2)

