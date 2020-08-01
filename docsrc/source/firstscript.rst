A first script
============================================================

Here is a first script that shows how DeerLab works. It simulates a noisy DEER data trace and then fits it.

We start by importing the required function/libraries:

.. code-block:: python

   from numpy import np
   from matplotlib.pyplot import plt
   from deerlab import *

We continue by generating a distance distribution consisting of a single Gaussian:

.. code-block:: python

   N = 201                            # number of points
   r = np.linspace(1.5,7,N)           # distance range, in nanometers
   P = dd_gauss(r,[3.5, 0.4])         # single-Gaussian distance distribution

Next, we calculate the background decay function due to a homogeneus 3D distribution of spins:

.. code-block:: python

   t = np.linspace(0,3,N)             # time axis, in microseconds
   conc = 100                         # spin concentration, in micromolar
   lam = 0.4                          # modulation depth
   B = bg_hom3d(t,conc,lam)           # homogeneous 3D background decay function

Next, we combine the distance distribution and the background into a full DEER signal and add some noise:

.. code-block:: python

   V = dipolarsignal(t,r,P,lam,B)     # DEER signal

   sig = 0.01                         # noise level
   Vexp = V + whitegaussnoise(V,sig)  # add noise

We can look at the result:

.. code-block:: python

   plt.plot(t,V,t,Vexp,t,(1-lam)*B)   # plotting
   plt.show()

Now that we have a noisy DEER trace, we fit it (in a single step) with a non-parametric distance distribution and a homogeneous 3D background.

.. code-block:: python

   Vfit,Pfit,Bfit,parfit,modeluq,stats = fitsignal(Vexp,t,r,'P',bg_hom3d)  # fitting

Finally, we plot the results

.. code-block:: python

   lamfit = parfit['ex']                         # fitted modulation depth
   
   # plotting
   plt.subplot(2,1,1)
   plt.plot(t,Vexp,t,Vfit,t,(1-lamfit)*Bfit)     # plot fitted model and background
   plt.subplot(2,1,2)
   plt.plot(r,P,r,Pfit)                          # plot fitted distribution
   plt.show()
