.. highlight:: matlab
.. _dipolarsignal:

*********************
:mod:`dipolarsignal`
*********************

Generation of dipolar signal from distance distribution

Syntax
=========================================

.. code-block:: matlab

    V = dipolarsignal(t,r,P)
    [V,S] = dipolarsignal(t,r,P)
    [V,S] = dipolarsignal(t,r,P,'Property',Value)


Parameters
    *   ``t`` - Time axis (N-array)
    *   ``r`` - Distance axis (M-array)
    *   ``P`` - Distance distribution (M-array)

Returns
    *   ``V`` - Form factor (N-array)
    *   ``S`` - Dipolar evolution function (N-array)


Description
=========================================

.. code-block:: matlab

    V = dipolarsignal(t,r,P)

Generates the noiseless form factor ``V`` on a time axis ``t`` from the distance distribution ``P`` in a distance axis ``r``.

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

       F = dipolarsignal(args,'ModDepth',0.4)


Background
    N-Array containing the background function of the form factor.

    *Default:* [*empty*]]

    *Example:*

    .. code-block:: matlab

       F = dipolarsignal(args,'Background',srtexp(t,param))

NoiseLevel
   Gaussian noise level (standard deviation)

    *Default:* ``0``

    *Example:*

    .. code-block:: matlab

        F = dipolarsignal(args,'NoiseLevel',0.05)

Overtones
    Array of RIDME overtone coefficients.

    *Default:* [*empty*]

    *Example:*

    .. code-block:: matlab

        F = dipolarsignal(args,'Overtones',[0.2 0.5 0.3])

Offset
    Vertical offset to multiply to the output signal

    *Default:* ``1``

    *Example:*

    .. code-block:: matlab

        F = dipolarsignal(args,'Offset', 1e3)


Examples
=========================================

.. code-block:: matlab


    t = linspace(-2,4,300);
    r = time2dist(t);
    P = onegaussian(r,[4 .3]);
    B = strexp(t,[0.15,3]);
    F = dipolarsignal(t,r,P,'NoiseLevel', 0.05,...
                            'ModDepth',0.4,...
                            'Background',B,...
                            'Offset',1000)

.. image:: ../images/dipolarsignal1.svg
