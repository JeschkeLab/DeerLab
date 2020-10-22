Quickstart tutorial
============================================================

This tutorial is intended as a overview of the fundamentals of DeerLab. If you are not familiar with DeerLab or 
dipolar EPR spectroscopy data analysis, then the :ref:`Beginner's Guide <beginners_guide>` might give you a better
insight on the basics for using DeerLab.  

This tutorial will assume that DeerLab is already installed on your computer. DeerLab always needs to be imported prior
to calling any of its functions as well as any other packages ::

   import deerlab as dl              # we always use 'dl' as the local name for DeerLab
   import numpy as np                # NumPy: vectors, matrices, linear algebra
   import matplotlib.pyplot as plt   # MatPlotLib: plotting

Loading experimental data
--------------------------




Simulating dipolar signals
---------------------------

DeerLab is provides the tools for simulating dipolar signals originating from different experiments. Note that these simulations 
are not based on spin dynamics simulations but rather on theoretical analytical treatments of such problems. 

Simulations in DeerLab can be easily achieved via the Fredholm integral which for a general dipolar signal reads

.. math:: V(t) = \int K(t,r)P(r) \mathrm{d}r

which in matrix form reads

.. math:: \boldsymbol{V} = \boldsymbol{K}\boldsymbol{P} 

To simulate a dipolar signal you need first to define a distance distribution. This distribution can be the result of a fit, come 
from molecular dynamics simulations, etc. or you can assume a certain shape for the distribution and use one of the parametric models 
available in DeerLab. For example, to simulate a Gaussian distance distribution use the ``dd_gauss`` model function

.. code-block:: python

   r = np.linspace(1.5,7,200)         # distance range, in nm
   rmean = 3.5 # nm
   sigma = 0.3 # nm
   P = dl.dd_gauss(r,[rmean, sigma])  # single-Gaussian distance distribution

Next, we calculate the background decay function due to a homogeneus 3D distribution of spins:

.. code-block:: python

   t = np.linspace(0,3,N)             # time axis, in microseconds
   conc = 100                         # spin concentration, in micromolar
   lam = 0.4                          # modulation depth
   B = dl.bg_hom3d(t,conc,lam)        # homogeneous 3D background decay function

Next, we combine the distance distribution and the background into a full 4-pulse DEER signal and add some noise:

.. code-block:: python

   K = dl.dipolarkernel(t,r,lam,B)       # DEER kernel
   V = K@P                               # DEER signal
   sig = 0.01                            # noise level
   Vexp = V + dl.whitegaussnoise(V,sig)  # add noise

We can look at the result:

.. code-block:: python

   plt.plot(t,V,t,Vexp)   # plotting
   plt.show()

Now that we have a noisy DEER trace, we fit itto it (in a single step) a model with a non-parametric distance distribution and a homogeneous 3D background.

.. code-block:: python

   fit = dl.fitsignal(Vexp,t,r,'P',dl.bg_hom3d,dl.ex_4pdeer)  # fitting

Finally, we plot the results:

.. code-block:: python

   # Plotting
   fit.plot()
