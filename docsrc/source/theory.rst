.. _theory:

Theory
=========================================

This section summarizes key aspects of the theory implemented in DeerLab.

Assumptions
-----------------------

The theory underlying DeerLab is not fully general and makes a series of assumptions about the sample and the experiment. The important ones are:

1. All spin labels are spin-1/2 or can be treated as such. High-spin systems are not handled by DeerLab.
2. There are no more than two spin labels on a protein (or other cluster). Multi-spin effects are not handled by DeerLab.
3. The unpaired electron spin density on the spin labels can be treated as localized in a single point. Delocalized spin systems are not handled by DeerLab.
4. All spins have essentially isotropic g values close to 2.00232. Spins with large g shifts and anisotropic g tensors are not handled by DeerLab.
5. There is no exchange coupling between any spins. The only interaction is through-space dipolar coupling. Exchange coupling is not handled by DeerLab.
6. The dipolar coupling between spins is in the weak-coupling regime, i.e. it is weaker than the difference between their resonance frequencies. Intermediate- and strong-coupling regimes are not handled by DeerLab.
7. All spins relax with the same phase memory time. Systems with rotamer-specific relaxation rates are not handled by DeerLab.
8. There is no orientation selection in the experiment. Orientation-selective DEER is not handled by DeerLab.

.. warning:: 
   If your sample and experiment do not satisfy all these assumptions, DeerLab will give incorrect results.


Distance distributions
-----------------------

In DeerLab, distance distributions are assumed to be normalized

.. math::
   \int_0^\infty P(r)\mr{d}r = 1 
   

Dipolar signals and kernels
-----------------------------------

A dipolar signal :math:`V(t)` is calculated from a distance distribution :math:`P(r)` with the Fredholm integral of the first kind:

.. math::

    V(t) = \int_0^\infty K(t,r)P(r)\mathrm{d}r

In this equation, :math:`K(t,r)` is the experiment-specific kernel representing how the signal is obtained from the distance distribution. Dipolar EPR experiments differ in the number and type of pulses, resulting in different number of dipolar pathways with different modulation amplitudes and refocusing times. In DeerLab, the signal from a general dipolar EPR experiment is modeled by the general kernel

.. math::
   K(t,r) = \left[\varLambda_0 + \sum_{p=1}^N \lambda_p K_0(n_p(t-T_p),r)\right]\cdot\prod_{p=1}^N B(n_p(t-T_p),\lambda_p)

Here, :math:`\varLambda_0` is the total amplitude of all unmodulated pathways, :math:`N` is the number of modulated pathways, :math:`\lambda_p` is the amplitude of the :math:`p`\ th modulated pathway, :math:`T_p` is the refocusing time of the :math:`p`\ th modulated pathway, and :math:`n_p` is the harmonic of the :math:`p`\ th pathway (most often, :math:`n=1`). `B` is the background decay function.

`K_0` is the elementary kernel for a single dipolar pathway with 100% modulation amplitude and unlimited excitation bandwidth. It is given by

.. math::

   K_0(t,r) =
   \int_0^1
   \cos\left((1-3\cos^2\theta) D r^{-3} t\right)
   \mathrm{d}\cos\theta

with the dipolar constant

.. math::

   D =
   \frac{\mu_0}{4\pi}
   \frac{(\mu_\mathrm{B}g_\mathrm{e})^2}{\hbar}
   \approx
   2\pi\cdot 52\,\mathrm{MHz\,nm^3}

The closed-form expression for :math:`K_0` is

.. math::

   K_0(t,r) = \frac{C(\xi)}{\xi}\cos(\omega_\perp t) + \frac{S(\xi)}{\xi} \sin(\omega_\perp t)

with :math:`\xi = \sqrt{6\omega_\perp t/\pi}` and the Fresnel cosine and sine integral functions

.. math::

   C(\xi) = \int_0^\xi \cos\left(\frac{\pi}{2}x^2\right)\mathrm{d}x
   \qquad
   S(\xi) = \int_0^\xi \sin\left(\frac{\pi}{2}x^2\right)\mathrm{d}x



For the common model used to analyze 4-pulse DEER data, the kernel is

.. math::
   K(t,r) = \left[(1-\lambda) + \lambda K_0(t,r)\right]\cdot B(t,\lambda)

This is a special case of the general kernel, with :math:`\varLambda_0 = 1-\lambda`, `N=1`, `\lambda_1=\lambda`, `T_1 = 0`, and `n_1=1`.


Discretization
-----------------------------

In DeerLab, all signals :math:`V` and distributions :math:`P` are represented as discretized vectors :math:`\vc{V}` and :math:`\vc{P}` over discretized time domain :math:`\vc{t}` and distance domain :math:`\vc{r}`. Their elements are

.. math::
   V_i = V(t_i)
   \qquad
   P_j = P(r_j)

Distance distribution vectors must have non-negative elements and are assumed to be normalized such that

.. math::
   \sum_j P_j \Delta r  = 1 

:math:`\Delta r` is the constant increment along the distance domain. DeerLab does not support non-linear distance vectors with non-constant increments.

All kernels :math:`K` are discretized accordingly to give kernel matrices :math:`\mx{K}` with elements

.. math::
   K_{ij} = K(t_i,r_j) \Delta r


With this, a signal is obtained from a distance distribution via

.. math::
   \vc{V} = \mx{K}\vc{P}



Least-squares fitting
-----------------------------

DeerLab uses dedicated least-squares solvers to fit models to data. The objective function and the solver depend on whether the distance distribution is parametric or parameter-free, and on whether there are background and experiment parameters to fit alongside the distance distribution.

Parametric distribution
.......................................

To fit a model with a parametric distance distribution to an experimental signal, DeerLab solves

.. math::

     \vc{\theta}_\mathrm{fit} =
     \argmin_{\vc{\theta}}
     \|\vc{V}_\mr{exp}-\mx{K}[\vc{\theta}]\vc{P}[\vc{\theta}]\|^2

where :math:`\vc{V}_\mr{exp}` indicates the experimental data and :math:`\vc{\theta}` is a vector of all parameters (distribution parameters, background parameters, experiment parameters). Various constrained least-squares solvers are implemented.

Parameter-free distribution
.......................................

To fit a model with a parameter-free distribution and no additional fitting parameters to an experimental signal, DeerLab implements several regularization approaches. The most common one is Tikhonov regularization. For this, the minimization problem is

.. math::

     \vc{P}_\mathrm{fit} =
     \argmin_{\vc{P}\ge0}
     \left(
     \|\vc{V}_\mr{exp}-\mx{K}\vc{P}\|^2
     +
     \alpha^2
     \|\mx{L}\vc{P}\|^2
     \right)

:math:`\alpha` is the regularization parameter, and :math:`\mx{L}` is the regularization operator matrix. DeerLab implements the linear non-negative least-squares solver FNNLS, as well as a few others. The :math:`\alpha` parameter can be optimized using a range of criteria, including L-curve, Akaike information criterion (AIC), and generalized cross validation (GCV).


To fit a  model with a parameter-free distance distribution and other parameters to an experimental signal, DeerLab solves

.. math::

    (\vc{\theta}_\mathrm{fit},\vc{P}_\mr{fit})
    =
    \argmin_{\vc{\theta},\vc{P}\ge0}
    \left(
    \|\vc{V}_\mr{exp}-\mx{K}[\vc{\theta}]\vc{P}\|^2
    +
    \alpha^2
    \|\mx{L}\vc{P}\|^2
    \right)

This problem is solved directly, i.e. both :math:`\vc{\theta}` ad :math:`\vc{P}` are fitted simultaneously. To achieve this, DeerLab implements a nested optimization approach that includes regularization.