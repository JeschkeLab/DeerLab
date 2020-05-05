A first script
============================================================

This section walks through a basic example DeerLab script that simulates a noisy DEER data trace and then fits it.

We start by generating a distance distribution consisting of a single Gaussian:

.. code-block:: matlab

   N = 201;                            % number of points
   r = linspace(1.5,7,N);              % distance range, in nanometers
   P = dd_gauss(r,[3.5 0.4]);          % single-Gaussian distance distribution

Next, we calculate the background decay function due to a homogeneus 3D distribution of spins:

.. code-block:: matlab

   t = linspace(0,3,N);                % time axis, in microseconds
   conc = 100;                         % spin concentration, in micromolar
   lambda = 0.4;                       % modulation depth
   B = bg_hom3d(t,conc,lambda);        % homogeneous 3D background decay function

Next, we combine the distance distribution and the background into a full DEER signal and add some noise:

.. code-block:: matlab

   V = dipolarsignal(t,r,P,lambda,B);  % DEER signal

   sig = 0.01;                         % noise level
   Vexp = V + whitegaussnoise(V,sig);  % add noise

We can look at the result:

.. code-block:: matlab

   plot(t,V,t,Vexp,t,(1-lambda)*B)     % plotting

Now that we have a noise DEER trace, we fit it (in a single step) with a parameter-free distance distribution and a homogeneous 3D background.

.. code-block:: matlab

   [Vfit,Pfit,Bfit,parfit] = fitsignal(Vexp,t,r,'P',@bg_hom3d);  % fitting

Finally, we plot the results

.. code-block:: matlab

   lambdafit = parfit.ex;                     % fitted modulation depth
   
   % plotting
   subplot(2,1,1)
   plot(t,Vexp,t,Vfit,t,(1-lambdafit)*Bfit);  % plot fitted model and background
   subplot(2,1,2)
   plot(r,P,r,Pfit);                          % plot fitted distribution

