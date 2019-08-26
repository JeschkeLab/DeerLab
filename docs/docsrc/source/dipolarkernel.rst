.. highlight:: matlab
.. _dipolarkernel:

*********************
:mod:`dipolarkernel`
*********************

Computes the dipolar interaction kernel for the linear transformation from distance-domain to time-domain.

"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
:mod:`K = dipolarkernel(t,r,B,...)`
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Parameters
    *   **t** - Time axis vector (N-array)
    *   **r** -  Distance axis vector (M-array)
    *   **B** -  Background function vector (N-array)
Returns
    *  **K** - Dipolar kernel (NxM-matrix)

Theory
=========================================
The dipolar kernel
:math:`K(t,r)`
is given by the following equation [1]_

.. math:: K(t,r)  =  \int_{0}^{\pi/2}\cos\left[(3\cos^2\theta -1)\frac{\mu_0\hbar\gamma_A\gamma_B}{4\pi r^3}t\right]\sin\theta d\theta

where
:math:`\omega_{dd}(r)`
is the dipolar modulation frequency,
:math:`\mu_0`
the permittivity of vacuum,
:math:`\hbar`
Planck's quantum of action and
:math:`\gamma_{A/B}`
are the gyromagnetic ratios of spin A and B, respectively.
This kernel describes the forward linear transformation from distance-domain to time-domain

.. math:: D(t) = \int_{0}^{\infty}K(t,r)P(r)dr.

For discretized data, this equation reads

.. math:: \mathbf{D} = \mathbf{K}\mathbf{P}\Delta\mathbf{r}.

therefore in order to account for the
:math:`\Delta\mathbf{r}`
factor, this can be introduced into the kernel to simplify the equation.

.. math:: \mathbf{K} = \mathbf{K}'\Delta\mathbf{r}.

It has been shown [2]_ that including the background function :math:`\mathbf{B}` into the kernel improves the solution of the inverse problem via least-squares methods. This is equivalent to defining the form factor directly in terms of the distance distribution and the kernel

.. math:: \mathbf{F} = (1-\lambda)  +  \lambda\mathbf{B}\mathbf{K}\mathbf{P}\Delta\mathbf{r}.

where :math:`\lambda` is the modulation depth. This again can be simplified by redifining the kernel as

.. math:: \mathbf{K} = (1-\lambda)  +  \lambda\mathbf{B}\mathbf{K}'\Delta\mathbf{r}

so that

.. math:: \mathbf{F} = \mathbf{K}\mathbf{P}.



Usage
=========================================
The function can be called as follows

.. code-block:: matlab

   K = dipolarkernel(t,r)
   K = dipolarkernel(t,r,[])

Computes ``K`` for the transformation to the dipolar evolution function from the time axis ``t`` and distance axis ``r``. The output kernel will not contain the background function.

.. Important:: Since  the dipolar kernel is normalized by
    :math:`\Delta r`
    , in order to obtain the correct time-domain signal, the distance distributions must be properly normalized by
    :math:`\Delta r`
    as well.

.. code-block:: matlab

    K = dipolarkernel(t,r,B)

If the background function ``B`` is specified, then it is included into the kernel function. Therefore, linear transformation of distance distributions using this kernel will yield the form factor with the background function.

.. Important:: If the dipolar kernel includes the background function, the form factor should not be background corrected.


Optional Arguments
=========================================
Optional arguments can be specified by parameter/value pairs. All property names are case insensitive and the property-value pairs can be passed in any order after the required input arguments have been passed..

.. code-block:: matlab

    K = dipolarkernel(args,'Property1',Value1,'Property2',Value2,...)

.. centered:: **Property Names & Descriptions**

KernelBType
    Specifies the way the background funcion ``B`` is introduced into the kernel:

    *   ``'full'`` - The input background is introduced unchanged

    *   ``'sqrt'`` - The square-root background is introduced

    *   ``'none'`` - background is not included in the kernel

    *Default:* ``'sqrt'``

    *Example:*

    .. code-block:: matlab

        K = dipolarkernel(args,'KernelBType','full') %Use background without changes

ExcitationBandwidth
    Excitation bandwith of the pulses in **MHz**. If specified, its value is used in the compensation of limited excitation bandwidth of the experimental pulses. If not specified infinite excitation bandwidth is assumed. The compensation for a given excitation bandwidth :math:`\Delta\omega` is taken into account by the approximation [3]_

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

KernelCalcMethod
    Specifies the way the kernel is computed numerically.


    *   ``'fresnel'`` - Employs Fresnel integrals for the kernel calculation (fast).

    *   ``'explicit'`` - Employs explicit powder averaging for the kernel calculation (slow).

    *Default:* ``'fresnel'``

    *Example:*

    .. code-block:: matlab

        K = dipolarkernel(args,'KernelCalcMethod','explicit')

Knots
    If the kernel is computed using the ``explicit`` powder averaging, this options specifies the number knots for the grid of powder orientations used for the powder averaging.

    *Default:* ``1001``

    *Example:*

    .. code-block:: matlab

        K = dipolarkernel(args,'KernelCalcMethod','explicit','Knots',2001)


References
=========================================

.. [1] Gunnar Jeschke, eMagRes, 2016, Vol 5: 1459–1476.
.. [2] Fábregas Ibáñez and Jeschke, to be published
.. [3] Banham et al., JMR 191, 2008, 202-218
