.. highlight:: matlab
.. _dipolarkernel_theory:

***********************
The Dipolar Kernel
***********************

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

For the case of the RIDME dipolar kernel, overtones up to the :math:`M^{th}` overtone can be included into the kernel as follows

    .. math:: K(t,r)  = \int_{0}^{\pi/2}\sum_{k=1}^M P_k\cos\left[(3\cos^2\theta -1)k\frac{\mu_0\hbar\gamma_A\gamma_B}{4\pi r^3}t\right]\sin\theta d\theta

    where :math:`P_k` are the overtone coefficients passed as arguments.

References
=========================================

.. [1] Gunnar Jeschke, eMagRes, 2016, Vol 5: 1459–1476.
.. [2] Fábregas Ibáñez and Jeschke, to be published
.. [3] Banham et al., JMR 191, 2008, 202-218
