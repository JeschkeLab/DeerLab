.. only:: html

    .. note::
        :class: sphx-glr-download-link-note

        Click :ref:`here <sphx_glr_download_auto_examples_plot_analyzing_pake_pattern.py>`     to download the full example code
    .. rst-class:: sphx-glr-example-title

    .. _sphx_glr_auto_examples_plot_analyzing_pake_pattern.py:


Analyzing the Pake pattern of a dipolar signal
============================================================================

A very basic example for displaying the Pake pattern of a given dipolar signal.


.. code-block:: python

    import numpy as np
    import matplotlib.pyplot as plt
    from deerlab import *








Generate a dipolar signal
-------------------------
Let's start by simulating a dipolar signal with some background and noise.

Prepare components


.. code-block:: python

    t = np.linspace(0,5,400)
    r = np.linspace(2,5,100)
    P = dd_gauss2(r,[3.5, 0.3, 0.2, 4, 0.2, 0.8])
    B = bg_exp(t,0.2)
    lam = 0.3
    K = dipolarkernel(t,r,lam,B)
    np.random.seed(0)
    V = K@P + whitegaussnoise(t,0.005)

    # Plot
    plt.plot(t,V,'k.')
    plt.tight_layout()
    plt.grid()
    plt.xlabel('Time [$\\mu s$]')
    plt.ylabel('V(t)')




.. image:: /auto_examples/images/sphx_glr_plot_analyzing_pake_pattern_001.png
    :alt: plot analyzing pake pattern
    :class: sphx-glr-single-img


.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none


    Text(42.597222222222214, 0.5, 'V(t)')



Prepare the signal
------------------
Since experimental dipolar signals contain the background, this must be fitted 
removed prior to Fourier transform.

First we proceed to fit the background function using some time-domain parametric 
model. In this example we will use an exponential function (``bg_exp``). 
Using the ``fitparamodel`` function we obtain the fitted background as well as 
the fitted modulation depth.


.. code-block:: python


    tstart = 3 # Time to start fitting background, in us
    mask = t>tstart
    # Model for the background component (1-lambda)*B
    def Bmodel(par):
        lam,kappa = par 
        B = (1 - lam)*bg_exp(t[mask],kappa)
        return B

    # Fit the background function
    fit = fitsignal(V,t,r,'P',bg_exp,ex_4pdeer,uqanalysis=False)
    Bfit = fit.B
    lam = fit.exparam
    kappa = fit.bgparam








Now we can use these fitted variables to isolate the dipolar evolution function 
from the primary data. Removal of the background via division leads to a noise 
increase at later times and thus to an approximation ``Vcorr`` of the real dipolar 
evolution function.


.. code-block:: python


    # "Correct" for the background and modulation depth
    Vcorr = (V/Bfit - (1 - lam))/lam

    plt.plot(t,Vcorr,'k.')
    plt.tight_layout()
    plt.grid(alpha=0.3)
    plt.xlabel('Time [$\\mu s$]')
    plt.ylabel('V(t)')




.. image:: /auto_examples/images/sphx_glr_plot_analyzing_pake_pattern_002.png
    :alt: plot analyzing pake pattern
    :class: sphx-glr-single-img


.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none


    Text(30.972222222222214, 0.5, 'V(t)')



Computing the Pake pattern
---------------------------

Now that the signal has the appropiate structure for Fourier transform it, 
we can call the ``fftspec`` function to obtained the Pake pattern.


.. code-block:: python


    # Compute spectrum
    nu,pake = fftspec(Vcorr,t,apodization=False)
 
     # %% [markdown]
    # In order to avoid truncation ripples in the Fourier spectrum and at the same 
    # time to compensate for the increase of noise, we recommend the use of apodization 
    # using the appropiate option in ``fftspec``.








Compute spectrum with apodization


.. code-block:: python

    nuapo,pakeapo = fftspec(Vcorr,t,apodization=False,mode='real')

    # Plot results
    plt.plot(nu,pake,'k',nuapo,pakeapo,'b',linewidth=1.5)
    plt.tight_layout()
    plt.grid(alpha=0.3)
    plt.xlim([-10, 10])
    plt.xlabel('Frequency [MHz]')
    plt.ylabel('Intensity [a.u.]')
    plt.legend(['Raw','Apodized'])




.. image:: /auto_examples/images/sphx_glr_plot_analyzing_pake_pattern_003.png
    :alt: plot analyzing pake pattern
    :class: sphx-glr-single-img


.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none


    <matplotlib.legend.Legend object at 0x000001ED98BDE198>



We do not need to worry about the zero-filling since ``fftspec`` takes care 
of setting it to twice the amount of points in the signal, to preserve all information. 
Adding more points will artificially increase the resolution of the Pake pattern.
The improvement will only be visual as no further information can be gained 
from additional zero-filling.


.. rst-class:: sphx-glr-timing

   **Total running time of the script:** ( 0 minutes  1.026 seconds)


.. _sphx_glr_download_auto_examples_plot_analyzing_pake_pattern.py:


.. only :: html

 .. container:: sphx-glr-footer
    :class: sphx-glr-footer-example



  .. container:: sphx-glr-download sphx-glr-download-python

     :download:`Download Python source code: plot_analyzing_pake_pattern.py <plot_analyzing_pake_pattern.py>`



  .. container:: sphx-glr-download sphx-glr-download-jupyter

     :download:`Download Jupyter notebook: plot_analyzing_pake_pattern.ipynb <plot_analyzing_pake_pattern.ipynb>`


.. only:: html

 .. rst-class:: sphx-glr-signature

    `Gallery generated by Sphinx-Gallery <https://sphinx-gallery.github.io>`_
