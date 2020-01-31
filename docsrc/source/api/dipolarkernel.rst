.. _dipolarkernel:

*********************
:mod:`dipolarkernel`
*********************

Computes the dipolar interaction kernel for the linear transformation from a distance-domain distribution to a time-domain dipolar signal.

-------------------------------


Syntax
=========================================

::

    K = dipolarkernel(t,r)
    K = dipolarkernel(t,r,lambda)
    K = dipolarkernel(t,r,lambda,B)
    K = dipolarkernel(t,r,'Property',Value)
    K = dipolarkernel(t,r,lambda,B,'Property',Value)


Parameters
    *   ``t``      - Time axis vector (*N*-element array)
    *   ``r``      -  Distance axis vector (*M*-element array)
    *   ``lambda`` - Modulation depth (scalar)
    *   ``B``      -  Background function vector (*N*-element array)
Returns
    *  ``K`` - Dipolar kernel (*NxM*-element matrix)

-------------------------------


Description
=========================================

.. code-block:: matlab

   K = dipolarkernel(t,r)

Computes ``K`` for the transformation to the dipolar evolution function from the time axis ``t`` and distance axis ``r``. This kernel will describe the transformation from the distance distribution `\mathbf{P}` to the dipolar evolution function `\mathbf{D}`, meaning


    .. math:: \mathbf{K}\mathbf{P}  = \mathbf{D}

Since the dipolar kernel is normalized by `\Delta r`, in order to obtain the correct time-domain signal via ``S=K*P``, the distance distribution must be properly normalized by `\Delta r` as well. By default, all distribution ``P`` returned by DeerAnalysis are already normalized.


-----------------------------


.. code-block:: matlab

    K = dipolarkernel(t,r,lambda)

If the modulation depth ``lambda`` is specified, then its information is included into the kernel function. This kernel will describe the transformation from the distance distribution to the form factor `\mathbf{F}` given by


    .. math:: \mathbf{K}\mathbf{P}  = \mathbf{F} = (1-\lambda) + \lambda \mathbf{D}


-----------------------------


.. code-block:: matlab

    K = dipolarkernel(t,r,lambda,B)

If the background ``B`` and modulation depth ``lambda`` variables are specified, then the background is included into the kernel function. This kernel will describe the transformation from the distance distribution to the primary signal `\mathbf{V}` given by

    .. math:: \mathbf{K}\mathbf{P}  = \mathbf{V} = [(1-\lambda) + \lambda \mathbf{D} ]\mathbf{B}


-------------------------------


Optional Arguments
=========================================


Optional arguments can be specified by parameter/value pairs. All property names are case insensitive and the property-value pairs can be passed in any order after the required input arguments have been passed..

.. code-block:: matlab

    K = dipolarkernel(r,t,'Property1',Value1,'Property2',Value2,...)
    K = dipolarkernel(r,t,lambda,'Property1',Value1,'Property2',Value2,...)
    K = dipolarkernel(r,t,lambda,B,'Property1',Value1,'Property2',Value2,...)

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

- ``'gValue'`` - Electron g-value
    Specifies the g-value of the electron spin center used to compute the dipolar frequencies from the given distance axis.

    *Default:* ``2.004602204236924``

    *Example:*

		.. code-block:: matlab

			K = dipolarkernel(args,'gValue',2.00) %Use experimental g-value

- ``'Method'`` - Numerical Kernel construction method
    Specifies the way the kernel is computed numerically.


    *   ``'fresnel'`` - Employs Fresnel integrals for the kernel calculation (fast).

    *   ``'explicit'`` - Employs explicit powder averaging for the kernel calculation (slow).

    *Default:* ``'fresnel'``

    *Example:*

		.. code-block:: matlab

			K = dipolarkernel(args,'Method','explicit')

- ``'Knots'`` - Number of powder orientations
    If the kernel is computed using the ``'explicit'`` powder averaging, this options specifies the number knots for the grid of powder orientations used for the powder averaging.

    *Default:* ``1001``

    *Example:*

    .. code-block:: matlab

        K = dipolarkernel(args,'Method','explicit','Knots',2001)

- ``'MultiPathway'`` - Dipolar multipathway kernel
    Parameters of the dipolar multi-pathway model. Passed as a cell array containing the following parameters:

                * ``etas``  - Time-shifts of the individual dipolar pathways. (N-vector)
                * ``lambdas``  - Amplitudes of the different pathways. The first element is the amplitude of the unmodulated component and the remaining elements are the amplitudes of the modulated components. (N+1-vector). 
                * ``Bmodel``  - Background model, must accept the time axis and the pathway amplitude as inputs. (function handle)

    *Default:* [*empty*]

    *Example:*

    .. code-block:: matlab

        lambdas = [0.5 0.3 0.2];
        etas = [0.0 5.4];
        Bmodel = @(t,lam) td_strexp(t,[kappa*lam d]);
        K = dipolarkernel(t,r,'MultiPathway',{lambdas etas Bmodel})
