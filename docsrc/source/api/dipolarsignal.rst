.. highlight:: matlab
.. _dipolarsignal:

*********************
:mod:`dipolarsignal`
*********************

Generate a dipolar signal from the distance distribution

-----------------------------



Syntax
=========================================

.. code-block:: matlab

    V = dipolarsignal(t,r)
    V = dipolarsignal(t,r,P)
    V = dipolarsignal(t,r,P,lam)
    V = dipolarsignal(t,r,P,lam,B)
    V = dipolarsignal(t,r,P,pathinfo)
    V = dipolarsignal(t,r,P,pathinfo,B)
    V = dipolarsignal(___,'Property',Value)


Parameters
    *   ``t`` - Time axis, in microseconds (*N*-element array)
    *   ``r`` - Distance axis, in nanometers (*M*-element array)
    *   ``P`` - Distance distribution (*M*-element array)
    *   ``lam`` - Modulation depth (scalar)
    *   ``pathinfo`` - Modulation pathway info array (Kx2 array)
    *   ``B`` - Backround (*N*-element array)

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


.. code-block:: matlab

    V = dipolarsignal(t,r,P,lam)

Include the modulation depth ``lam`` in the dipolar signal. If omitted, ``lam`` is set to 1.

-----------------------------

.. code-block:: matlab

    V = dipolarsignal(t,r,P,lam,B)

Include the background ``B`` in the dipolar signal. If omitted, ``B`` is set to all ``ones(size(t))``.

-----------------------------

.. code-block:: matlab

    V = dipolarsignal(t,r,P,pathinfo,B)

Compute the multi-pathway dipolar signal using the pathway specification in ``pathinfo``. See :ref:`dipolarkernel` for details.

-----------------------------


Additional Settings
=========================================

Additional settings can be specified via name-value pairs. All property names are case insensitive and the property-value pairs can be passed in any order after the required input arguments have been passed.

.. code-block:: matlab

    V = dipolarsignal(___,'Property1',Value1,'Property2',Value2,___)


- ``'NoiseLevel'`` - Level of noise on the signal
    Scalar value containing the desired standard deviation of a Gaussian noise vector 

    *Default:* ``0``

    *Example:*

		.. code-block:: matlab

			V = dipolarsignal(___,'NoiseLevel',0.05);

		.. Important::
			Each call of ``dipolarsignal`` will return a different noise realization. If you need a reproducible noise realization, seed MATLAB's random number generator with a specific integer seed ``k`` using ``rng(k)``.


- ``'Overtones'`` - RIDME overtone coefficients
    Array of RIDME overtone coefficients. The coefficients must be normalized, i.e. they must sum to unity.

    *Default:* [*empty*]

    *Example:*

		.. code-block:: matlab

			V = dipolarsignal(___,'Overtones',[0.2 0.5 0.3]);

- ``'g'`` - Electron g-value
    Specifies the g-value of the electron spin center used to compute the dipolar frequencies from the given distance axis.

    *Default:* ``2.00231930436256``

    *Example:*

		.. code-block:: matlab

			K = dipolarkernel(___,'g',2.005);   % Use experimental g-value

- ``'Scale'`` - Amplitude scale
    Vertical scale to multiply to the output signal

    *Default:* ``1``

    *Example:*

		.. code-block:: matlab

			V = dipolarsignal(___,'Scale', 1e3);

- ``'Phase'`` - IQ phase of the signal
    Scalar-valued phase of the complex-valued signal (in radians).

    *Default:* ``0``

    *Example:*

		.. code-block:: matlab

			V = dipolarsignal(___,'Phase', pi/2);
