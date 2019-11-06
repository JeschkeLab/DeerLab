.. highlight:: matlab
.. _dipolarsignal:

*********************
:mod:`dipolarsignal`
*********************

Generation of dipolar signals from distance distributions

Syntax
=========================================

.. code-block:: matlab

    V = dipolarsignal(t,r,P)
    [V,S] = dipolarsignal(t,r)
    [V,S] = dipolarsignal(t,r,P)
    [V,S] = dipolarsignal(t,r,P,'Property',Value)


Parameters
    *   ``t`` - Time axis (N-array)
    *   ``r`` - Distance axis (M-array)
    *   ``P`` - Distance distribution (M-array)

Returns
    *   ``V`` - Dipolar signal (N-array)
    *   ``S`` - Dipolar evolution function (N-array)


Description
=========================================

.. code-block:: matlab

    V = dipolarsignal(t,r,P)

Generates the noiseless form factor ``V`` on a time axis ``t`` from the distance distribution ``P`` in a distance axis ``r``.

.. code-block:: matlab

    V = dipolarsignal(t,r)

If no distribution ``P`` is provided, then the dipolar signal corresponding to a single distance ``r`` is computed.

Optional Arguments
=========================================
Optional arguments can be specified by parameter/value pairs. All property names are case insensitive and the property-value pairs can be passed in any order after the required input arguments have been passed..

.. code-block:: matlab

    V = dipolarsignal(args,'Property1',Value1,'Property2',Value2,...)


ModDepth
    Modulation depth of the form factor

    *Default:* ``1``

    *Example:*

    .. code-block:: matlab

       V = dipolarsignal(args,'ModDepth',0.4)


Background
    N-Array containing the background function of the form factor.

    *Default:* [*empty*]]

    *Example:*

    .. code-block:: matlab

       V = dipolarsignal(args,'Background',srtexp(t,param))

NoiseLevel
   Gaussian noise level (standard deviation)

    *Default:* ``0``

    *Example:*

    .. code-block:: matlab

        V = dipolarsignal(args,'NoiseLevel',0.05)

    .. Important::
       Each call of ``dipolarsignal`` will return a different noise realization. To set the output to a fixed noise realization, the random number generator must be fixed. In MATLAB this can be accomplished by calling ``rng(k)`` where ``k`` is some integer number.


Overtones
    Array of RIDME overtone coefficients.

    *Default:* [*empty*]

    *Example:*

    .. code-block:: matlab

        V = dipolarsignal(args,'Overtones',[0.2 0.5 0.3])

gValue
    Specifies the g-value of the electron spin center used to compute the dipolar frequencies from the given distance axis.

    *Default:* ``2.004602204236924``

    *Example:*

    .. code-block:: matlab

        K = dipolarkernel(args,'gValue',2.00) %Use experimental g-value

Scale
    Vertical scale to multiply to the output signal

    *Default:* ``1``

    *Example:*

    .. code-block:: matlab

        V = dipolarsignal(args,'Scale', 1e3)

Phase
    Phase of the complex-valued signal (in radians).

    *Default:* ``0``

    *Example:*

    .. code-block:: matlab

        V = dipolarsignal(args,'Phase', pi/2)

Interference
     Relative amplitude and time shift pairs of the dipolar interferences in multipulse-DEER experiments. The background model can be passed as a last argument to include the time-shifted backgrounds. 

    *Default:* [*empty*]

    *Example:*

    .. code-block:: matlab

        K = dipolarkernel(args,'Interference',[0.34 max(t)/2])
        K = dipolarkernel(args,'Interference',{0.34 max(t)/2 @td_strexp})


Examples
=========================================

.. code-block:: matlab


    t = linspace(-2,4,300);
    r = time2dist(t);
    P = onegaussian(r,[4 .3]);
    B = strexp(t,[0.15,3]);
    F = dipolarsignal(t,r,P,'NoiseLevel',0.05,...
                            'ModDepth',0.4,...
                            'Background',B,...
                            'Scale',1000)

.. image:: ../images/dipolarsignal1.svg
