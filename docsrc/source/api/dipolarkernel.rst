.. highlight:: matlab
.. _dipolarkernel:

*********************
:mod:`dipolarkernel`
*********************

Computes the dipolar interaction kernel `\mathbf{K}` for the linear transformation `\mathbf{S}=\mathbf{K}\mathbf{P}` from a distance-domain distribution `\mathbf{P}` to a time-domain signal `\mathbf{S}`.

Syntax
=========================================

.. code-block:: matlab

    K = dipolarkernel(t,r)
    K = dipolarkernel(t,r,lambda)
    K = dipolarkernel(t,r,lambda,B)
    K = dipolarkernel(t,r,'Property',Value)
    K = dipolarkernel(t,r,lambda,B,'Property',Value)


Parameters
    *   ``t`` - Time axis vector (N-array)
    *   ``r`` -  Distance axis vector (M-array)
    *   ``lambda`` - Modulation depth (scalar)
    *   ``B`` -  Background function vector (N-array)
Returns
    *  ``K`` - Dipolar kernel (NxM-matrix)

Description
=========================================

.. code-block:: matlab

   K = dipolarkernel(t,r)

Computes ``K`` for the transformation to the dipolar evolution function from the time axis ``t`` and distance axis ``r``. The output kernel will not contain the background function.

.. Note:: The function automatically detects the units of the inputs. The time axis ``t`` can be passed both in `\mu s` and `ns` and the distance axis ``r`` both in `nm` and Angstrom.

.. Important::
   Since the dipolar kernel is normalized by `\Delta r`, in order to obtain the correct time-domain signal with ``S=K*P``, the distance distribution must be properly normalized by `\Delta r` as well.

.. code-block:: matlab

    K = dipolarkernel(t,r,lambda,B)

If the background ``B`` and modulation depth ``lambda`` are specified, then the background is included into the kernel function. Therefore, linear transformation of distance distributions using this kernel will yield the form factor with the background function.

.. Important:: If the dipolar kernel includes the background function, the form factor should not be background corrected.


Optional Arguments
=========================================


Optional arguments can be specified by parameter/value pairs. All property names are case insensitive and the property-value pairs can be passed in any order after the required input arguments have been passed..

.. code-block:: matlab

    K = dipolarkernel(r,t,'Property1',Value1,'Property2',Value2,...)
    K = dipolarkernel(r,t,lambda,B,'Property1',Value1,'Property2',Value2,...)

ExcitationBandwidth
    Excitation bandwith of the pulses in **MHz**. If specified, its value is used in the compensation of limited excitation bandwidth of the experimental pulses. If not specified infinite excitation bandwidth is assumed. The compensation for a given excitation bandwidth :math:`\Delta\omega` is taken into account by the approximation [1]_

    .. math:: K(t,r,\Delta\omega)  = exp\left(-\frac{\omega_{dd}^2}{\Delta\omega^2}\right)K(t,r)

    *Default:* [*empty*]

    *Example:*

    .. code-block:: matlab

        K = dipolarkernel(args,'ExcitationBandwidth',50) %Correct for 50 MHz excitation bandwidth

OvertoneCoeffs
    1D-Array containing the overtone coefficients for RIDME experimens. If passed, the dipolar kernel overtones are calculated based on the passed coefficients. The coefficient values must be normalized. The kernel containing up to the :math:`K^{th}` overtone is constructed as follows

    .. math:: K(t,r)  = \int_{0}^{\pi/2}\sum_{k=1}^K P_k\cos\left[(3\cos^2\theta -1)k\frac{\mu_0\hbar\gamma_A\gamma_B}{4\pi r^3}t\right]\sin\theta d\theta

    where :math:`P_k` are the overtone coefficients passed as arguments.

    *Default:* [*empty*]

    *Example:*

    .. code-block:: matlab

        K = dipolarkernel(args,'OvertoneCoeffs',[0.4 0.2 0.4])

gValue
    Specifies the g-value of the electron spin center used to compute the dipolar frequencies from the given distance axis.

    *Default:* ``2.004602204236924``

    *Example:*

    .. code-block:: matlab

        K = dipolarkernel(args,'gValue',2.00) %Use experimental g-value

Method
    Specifies the way the kernel is computed numerically.


    *   ``'fresnel'`` - Employs Fresnel integrals for the kernel calculation (fast).

    *   ``'explicit'`` - Employs explicit powder averaging for the kernel calculation (slow).

    *Default:* ``'fresnel'``

    *Example:*

    .. code-block:: matlab

        K = dipolarkernel(args,'Method','explicit')

Knots
    If the kernel is computed using the ``explicit`` powder averaging, this options specifies the number knots for the grid of powder orientations used for the powder averaging.

    *Default:* ``1001``

    *Example:*

    .. code-block:: matlab

        K = dipolarkernel(args,'Method','explicit','Knots',2001)

FivePulseCoeff
    Two element array [A tshift] containing the relative amplitude of the 5-pulse DEER artefact and the time shift at which it appears. If not given, the time shift is set by default to half of tmax.

    *Default:* [*empty*]

    *Example:*

    .. code-block:: matlab

        K = dipolarkernel(args,'FivePulseCoeff',0.6) %5-pulse DEER artefact appearing at tau/2
        K = dipolarkernel(args,'FivePulseCoeff',[0.6 2]) %5-pulse DEER artefact appearing after 2us



References
=========================================

.. [1] Banham et al., JMR 191, 2008, 202-218
