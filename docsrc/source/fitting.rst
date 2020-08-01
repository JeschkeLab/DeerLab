Fitting
=========================================

DeerLab provides a wide range of functionality to analyze experimental dipolar EPR data using least-squares fitting. For fitting, the experimental data is assumed to be properly phased, shifted, and scaled.

The main fitting function is ``fitsignal``. It can fit either non-parametric and parametric distance distributions, and it fits the distribution, the background, and the modulation depth in a single step. It provides uncertainty estimates for all fitted quantities.

The older two-step workflow of fitting the background in a separate step is supported as well, but not recommended.


Picking the model
------------------------------------------

After preprocessing the experimental data, the next important step is to decide on the model to be fitted. This consists of four separate choices

(1) Choose a **distance range** ``r``. The distance range is an important choice, as any distance distibution is truncated to this range. The lower limit of the distance range is determined by the bandwidth of the pulses, and also on the time increment. Typically, 1.5 nm is reasonable. The upper limit depends on the length of the experimental time trace. The number of points in ``r`` is usually set equal to the number of time points.

(2) Choose a **distribution model**. Generally, a non-parametric distribution is preferred. If there a reasons to believe the distances distribution has a specific shape, use the associated :ref:`parametric distance distibution model<modelsref_dd>`.

(3) Choose a **background model**. Typically, a homogeneous 3D background (``bg_hom3d``) is sufficient, but other :ref:`background models<modelsref_bg>` could be needed.

(4) Select an **experiment model**. If analyzing a standard 4-pulse DEER trace without 2+1 component at the end, use ``ex_4pdeer``. If the 2+1 components should be fitted as well, use ``ex_ovl4pdeer``. There are :ref:`experimental models<modelsref_ex>` for more complicated signals.


Non-parametric distribution
------------------------------------------

To fit a signal with a non-parametric distance distribution, the simplest possible call is

.. code-block:: python

    Vfit,Pfit,Bfit,parfit,paruq,moduq,stats = fitsignal(Vexp,t,r)

This fits the signal in ``Vexp`` (defined over the time axis ``t``) with a non-parametric distance distribution defined over distance vector ``r``, a homogeneous 3D spin distribution for the background, and a modulation depth (assuming a standard 4-pulse DEER).

If ``fitsignal`` is not asked for an output, the results of the fit, including uncertainty, are plotted in a figure and printed to the command window.

To obtain the fitted signal, use

.. code-block:: python

    Vfit,Pfit,Bfit,parfit,paruq,moduq,stats = fitsignal(Vexp,t,r)

With additional inputs, you can specify which :ref:`distance distribution model <modelsref_dd>`, :ref:`background model <modelsref_bg>`, and :ref:`experiment model <modelsref_ex>` to use for fitting:

.. code-block:: python

    Vfit,Pfit,Bfit,parfit,paruq,moduq,stats = fitsignal(Vexp,t,r,'P',bg_hom3d,ex_4pdeer)

Here, ``'P'`` specifies that a non-parametric distribution (defined over ``r``) should be fitted. ``bg_hom3d`` is a function handle to the background model function. If you want to fit a signal that doesn't have a background, use  ``'none'`` instead.  ``ex_4pdeer`` is a function handle to the particular experiment model. ``ex_4pdeer`` consists of the main dipolar modulation function centered at time zero with an ampltitude ``lambda`` plus a constant offset of amplitude ``1-lambda``. If you want to fit a simple dipolar evolution function with 100% modulation depth, use ``'none'`` as the experiment model.

``fitsignal`` uses a least-squares fitting algorithm to determine the optimal distribution, background parameters, and experiment parameters that fit the experiment data. To determine the non-parametric distribution, it internally uses Tikhnonov regularization with a regularization parameter optimized using the Akaike Information Criterion (AIC). These settings can be changed:

.. code-block:: python

   regtype = 'tv'     # use total variation instead of Tikhonov regularization
   alpha = 0.8        # manually set regularization parameter
   Vfit,Pfit,Bfit,parfit,paruq,moduq,stats = fitsignal(Vexp,t,r,'P',bg_hom3d,ex_4pdeer,regtype=regtype,regparam=alpha)

``fitsignal`` can return additional outputs: the fitted background (``Bfit``), the fitted distribution (``Pfit``), the fitted parameters (``parfit``), and confidence intervals for all parameters (``parci``):

.. code-block:: python

    Vfit,Pfit,Bfit,parfit,paruq,moduq,stats = fitsignal(Vexp,t,r,'P',bg_hom3d,ex_4pdeer)

Parametric distributions
----------------------------------

To use a parametric distance distribution model, provide ``fitsignal`` with a function handle to the :ref:`distance distribution models<modelsref_dd>` instead of ``'P'``. For example:

.. code-block:: python

    Vfit,Pfit,Bfit,parfit,paruq,moduq,stats = fitsignal(Vexp,t,r,dd_gauss2,bg_hom3d,ex_4pdeer)

This will fit a two-Gauss distribution over ``r``. The fitted distribution parameters are returned in ``parfit``, and the corresponding distribution in ``Pfit``.