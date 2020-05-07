Fitting
=========================================

DeerLab provides a wide range of functionality to analyze experimental dipolar EPR data using least-squares fitting. For fitting, the experimental data is assumed to be properly phased, shifted, and scaled.

The main fitting function is ``fitsignal``. It can fit either parameter-free and parametric distance distributions, and it fits the distribution, the background, and the modulation depth in a single step. It provides uncertainty estimates for all fitted quantities.

The older two-step workflow of fitting the background in a separate step is supported as well, but not recommended.


Picking the model
------------------------------------------

After preprocessing the experimental data, the next important step is to decide on the model to be fitted. This consists of four separate choices

(1) Choose a **distance range** ``r``. The distance range is an important choice, as any distance distibution is truncated to this range. The lower limit of the distance range is determined by the bandwidth of the pulses, and also on the time increment. Typically, 1.5 nm is reasonable. The upper limit depends on the length of the experimental time trace. One possible heuristic is ``rmax = 6*(tmax/2)^(1/3)``. The number of points in ``r`` is usually set equal to the number of time points.

(2) Choose a **distribution model**. Generally, a parameter-free distribution is preferred. If there a reasons to believe the distances distribution has a specific shape, use the associated :ref:`parametric distance distibution model<modelsref_dd>`.

(3) Choose a **background model**. Typically, a homogeneous 3D background (``bg_hom3d``) is sufficient, but other :ref:`background models<modelsref_bg>` could be needed.

(4) Select an **experiment model**. If analyzing a standard 4-pulse DEER trace without 2+1 component at the end, use ``ex_4pdeer``. If the 2+1 components should be fitted as well, use ``ex_ovl4pdeer``. There are :ref:`experimental models<modelsref_ex>` for more complicated signals.


Parameter-free distribution
------------------------------------------

To fit a signal with a parameter-free distance distribution, the simplest possible call is

.. code-block:: matlab

    fitsignal(Vexp,t,r)

This fits the signal in ``Vexp`` (defined over the time axis ``t``) with a parameter-free distance distribution defined over distance vector ``r``, a homogeneous 3D spin distribution for the background, and a modulation depth (assuming a standard 4-pulse DEER).

If ``fitsignal`` is not asked for an output, the results of the fit, including uncertainty, are plotted in a figure and printed to the command window.

To obtain the fitted signal, use

.. code-block:: matlab

    Vfit = fitsignal(Vexp,t,r)

With additional inputs, you can specify which :ref:`distance distribution model <modelsref_dd>`, :ref:`background model <modelsref_bg>`, and :ref:`experiment model <modelsref_ex>` to use for fitting:

.. code-block:: matlab

    Vfit = fitsignal(Vexp,t,r,'P',@bg_hom3d,@ex_4pdeer)

Here, ``'P'`` specifies that a parameter-free distribution (defined over ``r``) should be fitted. ``@bg_hom3d`` is a function handle to the background model function. If you want to fit a signal that doesn't have a background, use  ``'none'`` instead.  ``@ex_4pdeer`` is a function handle to the particular experiment model. ``ex_4pdeer`` consists of the main dipolar modulation function centered at time zero with an ampltitude ``lambda`` plus a constant offset of amplitude ``1-lambda``. If you want to fit a simple dipolar evolution function with 100% modulation depth, use ``'none'`` as the experiment model.

``fitsignal`` uses a least-squares fitting algorithm to determine the optimal distribution, background parameters, and experiment parameters that fit the experiment data. To determine the parameter-free distribution, it internally uses Tikhnonov regularization with a regularization parameter optimized using the Akaike Information Criterion (AIC). These settings can be changed:

.. code-block:: matlab

   regtype = 'tv';     % use total variation instead of Tikhonov regularization
   alpha = 0.8;        % manually set regularization parameter
   Vfit = fitsignal(Vexp,t,r,'P',@bg_hom3d,@ex_4pdeer,'RegType',regtype,'RegParam',alpha)

``fitsignal`` can return additional outputs: the fitted background (``Bfit``), the fitted distribution (``Pfit``), the fitted parameters (``parfit``), and confidence intervals for all parameters (``parci``):

.. code-block:: matlab

    [Vfit,Bfit,Pfit,parfit,parci] = fitsignal(Vexp,t,r,'P',@bg_hom3d,@ex_4pdeer)

Parametric distributions
----------------------------------

To use a parametric distance distribution model, provide ``fitsignal`` with a function handle to the :ref:`distance distribution models<modelsref_dd>` instead of ``'P'``. For example:

.. code-block:: matlab

    [Vfit,Bfit,Pfit,parfit,parci] = fitsignal(Vexp,t,r,@dd_gauss2,@bg_hom3d,@ex_4pdeer)

This will fit a two-Gauss distribution over ``r``. The fitted distribution parameters are returned in ``parfit``, and the corresponding distribution in ``Pfit``.


Background fitting
------------------------------------------

Although a fitted background can be obtained directly in a one-step analysis from ``fitsignal`` (see above), for compatibility with other programs such as DeerAnalysis, DeerLab also allows to fit a background function separately to a dipolar signal. This is the first step in the traditional two-step analysis.

To fit a homogeneous 3D background and a modulation amplitude to a dipolar signal, use

.. code-block:: matlab

   [Bfit,lamfit] = fitbackground(t,V,@bg_hom3d)

Other background models can be fitted as well. If the modulation depth is known, it can be provided to ``fitbackground``. 
``fitbackground`` fits the background to the tail of the dipolar signal. The fitting range can be adjusted by providing a starting and an end point:

.. code-block:: matlab

   [Bfit,lamfit] = fitbackground(t,V,@bg_hom3d,[tstart tend]);

If these are not provided, ``fitbackground`` tries to determine an optimal range automatically.

Once a background and a modulation ampltiude are obtained, the second step in the traditional two-step approach is either to divide out the background and the modulation amplitude and then fit the remaining signal ``D``:

.. code-block:: matlab

   D = (V./Bfit - (1-lamfit))/lamfit;
   K0 = dipolarkernel(t,r);
   Pfit = fitregmodel(D,K0,r);

or to include the background decay and modulation amplitude into the kernel:

.. code-block:: matlab

   K = dipolarkernel(t,r,lamfit,Bfit);
   Pfit = fitregmodel(V,K,r);

For long signals with many oscillations or a long stretch of background-only decay, the two-step analysis and the one-step analysis give similar results. For truncated signals that are only a few oscillation periods long (or less), the one-step analysis with ``fitsignal`` is significantly more robust.
