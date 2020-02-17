.. _dipolarkernel:

*********************
:mod:`dipolarkernel`
*********************

Computes the dipolar interaction kernel matrix ``K`` that transforms a distance distribution ``P`` to a time-domain dipolar signal ``D`` via ``D = K*P``.

-------------------------------


Syntax
=========================================

::

    K = dipolarkernel(t,r)
    K = dipolarkernel(t,r,lambda)
    K = dipolarkernel(t,r,lambdaT0)
    K = dipolarkernel(t,r,lambda,B)
    K = dipolarkernel(t,r,lambdaT0,B)
    K = dipolarkernel(...,'Property',Value)


Parameters
    *   ``t``        - Time axis vector (*N*-element array)
    *   ``r``        - Distance axis vector (*M*-element array)
    *   ``lambda``   - Modulation depth (scalar)
    *   ``lambdaT0`` - Array of modulation depths and refocusing times (*mx2* array)
    *   ``B``        - Background, either vector of values (*N*-element array) or function handle
Returns
    *  ``K`` - Dipolar kernel (*NxM* array)

-------------------------------


Description
=========================================

.. code-block:: matlab

   K = dipolarkernel(t,r)

Computes the kernel matrix ``K`` for the time axis ``t`` and distance axis ``r``. ``K`` will describe the transformation from the distance distribution `\mathbf{P}` to the dipolar evolution function `\mathbf{D}`

    .. math:: \mathbf{K}\mathbf{P}  = \mathbf{D}

Since the dipolar kernel is normalized by `\Delta r`, in order to obtain the correct time-domain signal via ``D=K*P``, the distance distribution must be properly normalized by `\Delta r` as well. By default, all distributions ``P`` returned by DeerLab model distribution functions are already normalized.


-----------------------------


.. code-block:: matlab

    K = dipolarkernel(t,r,lambda)

If the modulation depth ``lambda`` is specified, then its information is included into the kernel function. This kernel will describe the transformation from the distance distribution to the form factor `\mathbf{F}` given by


    .. math:: \mathbf{K}\mathbf{P}  = \mathbf{F} = (1-\lambda) + \lambda \mathbf{D}


-----------------------------


.. code-block:: matlab

    K = dipolarkernel(t,r,lambda,B)

If the background ``B`` and modulation depth ``lambda`` variables are specified, then the background is included into the kernel function. ``B`` can be either an array with the precalculated background decay, or a function handle. ``K`` will describe the transformation from the distance distribution to the primary signal `\mathbf{V}` given by

    .. math:: \mathbf{K}\mathbf{P}  = \mathbf{V} = [(1-\lambda) + \lambda \mathbf{D} ]\mathbf{B}


-------------------------------


.. code-block:: matlab

    K = dipolarkernel(t,r,lambdaT0)
    K = dipolarkernel(t,r,lambdaT0,B)

For a multi-pathway DEER signal (e.g, 4-pulse DEER with 2+1 contribution; 5-pulse DEER with 4-pulse DEER residual signal), ``lambdaT0`` contains a list of modulation depths (amplitudes) and refocusing times (in microsecods) for all modulated pathway signals. Each row of ``lambdaT0`` needs two values: one amplitude and one time. The unmodulated amplitudes is calculated such that the sum over all amplitudes equals 1.

.. code-block:: matlab

    lambda = [0.7 0.2];
    T0 = [0 4];
    lambdaT0 = [lambda(:) T0(:)];
    K = dipolarsignal(t,r,lambdaT0);



Additional Settings
=========================================


Additional settings can be specified by parameter/value pairs. All property names are case insensitive and the property-value pairs can be passed in any order after the required input arguments have been passed..

.. code-block:: matlab

    K = dipolarkernel(...,'Property1',Value1,'Property2',Value2,...)

- ``'ExcitationBandwidth'`` - Excitation bandwith of the pulses in **MHz**. 
    If specified, its value is used in the compensation of limited excitation bandwidth of the experimental pulses. If not specified infinite excitation bandwidth is assumed. The compensation for a given excitation bandwidth :math:`\Delta\omega` is taken into account by the approximation

    .. math:: K(t,r,\Delta\omega)  = exp\left(-\frac{\omega_{dd}^2}{\Delta\omega^2}\right)K(t,r)

    *Default:* [*empty*]

    *Example:*

    .. code-block:: matlab

        K = dipolarkernel(args,'ExcitationBandwidth',50) %Correct for 50 MHz excitation bandwidth

- ``'OvertoneCoeffs'`` - RIDME overtone coefficients
    1D-Array containing the overtone coefficients for RIDME experimens. If passed, the dipolar kernel overtones are calculated based on the passed coefficients. The coefficient values must be normalized. The kernel containing up to the :math:`K^{th}` overtone is constructed as follows

    .. math:: K(t,r)  = \int_{0}^{\pi/2}\sum_{k=1}^K P_k\cos\left[(3\cos^2\theta -1)k\frac{\mu_0\hbar\gamma_A\gamma_B}{4\pi r^3}t\right]\sin\theta d\theta

    where :math:`P_k` are the overtone coefficients passed as arguments.

    *Default:* [*empty*]

    *Example:*

		.. code-block:: matlab

			K = dipolarkernel(args,'OvertoneCoeffs',[0.4 0.2 0.4])

- ``'g'`` - Electron g-value
    Specifies the g-value of the electron spin center used to compute the dipolar frequencies from the given distance axis.

    *Default:* free-electron g value

    *Example:*

		.. code-block:: matlab

			K = dipolarkernel(args,'g',2.01)

- ``'Method'`` - Numerical Kernel construction method
    Specifies the way the kernel is computed numerically.


    *   ``'fresnel'`` - Uses Fresnel integrals for the kernel calculation (fast).

    *   ``'grid'`` - Uses powder averaging over a grid of orientations for the kernel calculation (slow).

    *Default:* ``'fresnel'``

    *Example:*

		.. code-block:: matlab

			K = dipolarkernel(args,'Method','grid')

- ``'nKnots'`` - Number of powder orientations
    If the kernel is computed using ``'grid'``, this options specifies the number knots for the grid of powder orientations used for the powder averaging.

    *Default:* ``1001``

    *Example:*

    .. code-block:: matlab

        K = dipolarkernel(args,'Method','grid','nKnots',2001)
