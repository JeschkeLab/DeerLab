Fitting
=========================================

DeerLab provides a wide range of functionality to analyze experimental dipolar EPR data. The main fitting function is ``fitsignal``.  The experimental data is assumed to be properly phased, shifted, and scaled. ``fitsignal`` can fit either parameter-free and parametric distance distributions, and it fits the distribution, the background, and the modulation depth in a single step. It provides uncertainty estimates for all fitted quantities.

The older two-step workflow of fitting the background in a separate step is supported as well, but not recommended.


Model selection
------------------------------------------

After preprocessing the data, the next important step is to decide on the model to be fitted. This consists of four steps

1. Choose a distance range ``r``. If unclear, use ``r = time2dist(t)``, which will estimate an accessible distance range based on the time axis. The distance range is an important choice, as any distance distibution is truncated to this range.
2. Choose a distribution model, either parameter-free or parameterized. Generally, a parameter-free model is preferred.
3. Choose a background model. Typically, a homogeneous 3D background (``bg_hom3d``) is sufficient, but other background models could be needed.
4. Select an experiment model. If analyzing a standard 4-pulse DEER trace without 2+1 component at the end, use ``ex_4pdeer``.


Parameter-free distribution
------------------------------------------

The simplest possible call is

.. code-block:: matlab

    fitsignal(Vexp,t,r)

This fits the signal in ``Vexp`` (defined over the time axis ``t``) with a parameter-free distance distribution defined over distance vector ``r``, a homogeneous 3D spin distribution for the background, and a modulation depth (assuming a standard 4-pulse DEER). The results of the fit, including uncertainty, are plotted.

To obtain the fitted signal, use

.. code-block:: matlab

    Vfit = fitsignal(Vexp,t,r)

With additional inputs, you can specify which :ref:`distance distribution model <ddmodels>`, :ref:`background model <bgmodels>`, and :ref:`experiment model <exmodels>` to use for fitting:

.. code-block:: matlab

    Vfit = fitsignal(Vexp,t,r,'P',@bg_hom3d,@ex_4pdeer)

Here, ``'P'`` specifies that a parameter-free distribution (defined over ``r``) should be fitted.

``@bg_hom3d`` is a function handle to the :ref:`background model <bgmodels>` to be used. If you want to fit a signal that doesn't have a background, use  ``'none'`` instead. 

``@ex_4pdeer`` is a function handle to the particular :ref:`experiment model <exmodels>`. ``ex_4pdeer`` consists of the main dipolar modulation function centered at time zero with an ampltitude ``lambda`` plus a constant offset of amplitude ``1-lambda``.

``fitsignal`` can return additional outputs such as the fitted background(``Bfit``), the fitted distribution (``Pfit``), the fitted parameters (``parfit``), and confidence intervals for the parameters (``parci``):

.. code-block:: matlab

    [Vfit,Bfit,Pfit,parfit,parci] = fitsignal(Vexp,t,r,'P',@bg_hom3d,@ex_4pdeer)

Parameterized distributions
----------------------------------

If you prefer to use a parameterized distance distribution model, provide a function handle to one of the many distance distribution models instead. For example:

.. code-block:: matlab

    [Vfit,Bfit,Pfit,parfit,parci] = fitsignal(Vexp,t,r,'P',@bg_hom3d,@ex_4pdeer)

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
   K = dipolarkernel(t,r);
   Pfit = fitregmodel(D,K,r);

or to include the background and modulation amplitude into the kernel:

.. code-block:: matlab

   K = dipolarkernel(t,r,lamfit,Bfit);
   Pfit = fitregmodel(V,K,r);

For long signals with many oscillations or a long stetch of background-only decay, the two-step analysis lead to similar results as the one-step analysis. For truncated signals that are only a few oscillation periods long (or less), one-step analysis is significantly more robust.
