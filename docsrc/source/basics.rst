Basics
=========================================

DeerLab relies on a few central quantities. They include distance distributions and model functions, background and their model functions, time-domain signals, and kernels. They are described in this section.

All functions in DeerLab use the same units: All distances are in units of **nanometers**, and all times in units of **microseconds**.

Distance distributions
------------------------------------------

A distance distribution between two spins is represented by a pair of vectors: a distance vector ``r`` (in nanometers) and a vector of densities ``P`` (in inverse nanometers). ``r`` can be a linearly or non-linearly increasing vector, and it cannot have negative values. The elements in ``P`` must be nonnegative. ``P[i]`` is the distribution value at ``r[i]``. Outside the range of ``r``, ``P`` is assumed to be zero, i.e. the distribution is truncated to the range ``r``. ``P`` is normalized such that the integral over the range of the provided ``r`` equals one:

.. code-block:: python

        np.trapz(P,r)                 # integral of P over r, gives 1


DeerLab distinguishes between **non-parametric** and **parametric** distance distributions.

**Non-parametric distance distributions** are the more general from. They have no particular shape and are represented by vectors ``P`` and ``r``. For example, you can generate ``P`` and ``r`` by an external program for spin label rotamer modeling. Non-parametric distance distributions are returned in least-squares fitting. In least-squares fitting, non-parametric distance distributions are preferred over parametric distance distributions, since they make fewer assumptions about the distribution and are more flexible. They introduce less bias.

**Parametric distance distributions** have specific shapes that depend on a few parameters. DeerLab provides many parametric :ref:`distance distribution model functions <modelsref_dd>`. All these functions start with the prefix ``dd_`` (``dd`` stands for "distance distibution"). They take a distance vector ``r`` and a parameter vector ``param`` as inputs and return the distance distribution as a vector ``P``. Here is an example:

.. code-block:: python
   
        r = np.linspace(1.5,6,201)   # distance, in nanometers
        P = dd_gauss(r,[3 0.5])      # Gaussian distribution, center 3 nm, fwhm 0.5 nm
        plt.plot(r,P)

The first line generates a distance axis ``r`` as a row vector. The second line generates the distance distribution ``P`` as a vector defined over ``r``. The function ``dd_gauss`` is a distance distribution model function provided by DeerLab. It takes as inputs the distance axis and a two-element parameter array (center and full width at half maximum) and returns the calculated distance distribution consisting of one Gaussian, truncated to the range of ``r`` and normalized.

To programmatically get information on a particular :ref:`distance distribution model functions <modelsref_dd>` and its parameters, call the function without input arguments

.. code-block:: python

        info = dd_gauss()        # obtain information on model and parameters

DeerLab provides a wide range of :ref:`parametric distribution models<modelsref_dd>` that fall into several groups. There are models that are based on basis functions and linear combinations thereof. Gaussian and 3D-Rice functions are the most important ones. To use a single Gaussian distribution, use :ref:`dd_gauss`. To combine 2, 3, etc. Gaussians distributions, use :ref:`dd_gauss2`, :ref:`dd_gauss3`, etc. Similarly, for 3D-Rice distributions, use :ref:`dd_rice`, :ref:`dd_rice2`, etc. A second group of distance distribution models represent three-dimensional disordered segmented objects such as proteins and other polymers: :ref:`dd_wormchain`, :ref:`dd_wormgauss`, :ref:`dd_randcoil`. A third group provides models for distributions of labels in simple confined spaces such as spheres and spherical shells. Finally, DeerLab provides a series of distance distribution models that represent simple geometric shapes. These models are worthless for practical purposes, but are useful as toy distributions for methods development: :ref:`dd_circle`, :ref:`dd_triangle`, and :ref:`dd_uniform`.


.. _bgmodels:


Backgrounds
--------------------------------------------

DeerLab includes a collection of parametric models that can be used to model the background or inter-molecular signal, i.e. the signal due to randomly distributed spins in the sample that are not part of the spin-labeled protein or object. All these :ref:`background model functions <modelsref_bg>` start with the prefix ``bg_`` and have the same calling syntax. The inputs a time axis vector ``t`` (in microseconds), a parameter vector ``param``, and a modulation amplitude ``lambda`` (between 0 and 1) . The length of ``param``, and the meaning of the elements, depends on the particular model. If ``lambda`` is not provided, it is set to one. The output is a background decay vector ``B``, defined over ``t``.

Here is an example:

.. code-block:: python

        t = np.linspace(-0.1,3,201) # time, in microseconds
        lam = 0.4                   # modulation depth
        conc = 200                  # spin concentration, in uM
        B = bg_hom3d(t,conc,lam)    # homogeneous 3D background
        plt.plot(t,B)

The first line generate the desired time axis. The second line gives the modulation depth, and the third gives the spin concentration (in micromolar). Both are inputs to the background function ``bg_hom3d``, which calculates a decay due to a homogeneous three-dimensional distribution of spins and returns it in ``B``. 

To get information on the model and its parameters, call the function without inputs:

.. code-block:: python

        info = bg_hom3d()        # obtain information on model and parameters


DeerLab's :ref:`background models<modelsref_bg>` fall into two categories, physical and phenomenological. **Physical models** describe particular distributions of spin labels in space. These models depend on physical parameters such as spin concentration, exclusion distances, and dimensionality. The most common one is :ref:`bg_hom3d`, which describes the signal due to a homogeneous three-dimensional distribution of spins of a given concentration. A homogeneous distribution in a fractal dimensions is available with :ref:`bg_homfractal`, and excluded-volume effects can be modelled using :ref:`bg_hom3dex`. **Phenomenological models** represent various mathematical functions that are intended to mimick the background decay, without reference to a particular spatial distribution of spins. The parameters of these models do no have direct physical meaning. In general, it is preferable to use the physical instead of phenomenological models.


.. _exmodels:

Experiments
-------------------------------------------------------------

DeerLab supports a wide range of dipolar EPR experiments. Experiments differ in the number of dipolar modulation components and their refocusing times. For each type of supported dipolar EPR experiment, there is a dedicated :ref:`experiment model function<modelsref_ex>` starting with ``ex_``. These functions take as inputs the time axis ``t`` and an array of parameters characterizing the experiment. As output, they return an array containing information about the dipolar pathways of the experiment model.

For example, the model function representing the typical model for a 4-pulse DEER signal is ``ex_4pdeer``:

.. code-block:: python

        t = np.linspace(0,3,151)
        lam = 0.3;
        pathways = ex_4pdeer(t,lam)

The returned output is

.. code-block:: python

    pathways =[[0.7,  NaN],
               [0.3,   0 ]]

Each row of this array holds information about one pathway. The first column is modulation amplitude, and the second column is the refocusing point. In the above example, the first row shows a pathway with amplitude 0.7 and no refocusing time, indicating that it represents the unmodulated contribution. The pathway of the second row shows amplitude of 0.3 and refocusing time 0, i.e. this is the primary dipolar pathway.

Kernel matrices
--------------------------------------------

One of the core functions of DeerLab is ``dipolarkernel``. It provides the kernel that provides the connection between the distance distribution and the time-domain signal.

.. code-block:: python

    t = np.linspace(0,3,301)     # time axis, in us
    r = np.linspace(2,7,301)     # distance axis, in nm
    K0 = dipolarkernel(t,r)      # dipolar kernel matrix

To obtain the time-domain signal due to a distribution ``P``, use

.. code-block:: python
    
    V = K0@P

The above is the most elementary kernel, giving a time-domain signal without any background decay, and with a single dipolar evolution function centered at time zero and with modulation depth of 1.

The kernel can also include the background and the modulation depth. Then, the multiplication of ``P`` by ``K`` will return the complete time-domain signal. Here is an example:

.. code-block:: python

    lam = 0.4
    B = bg_hom3d(t,200,lam)
    K = dipolarkernel(t,r,lam,B)
    V = K@P

The function ``dipolarkernel`` also has options to add an excitation bandwidth limitation, to select the internal calculation method, and more.

It is not necessary to precompute the background decay. Instead, provide ``dipolarkernel`` with a lambda-callable to a function that takes only time and modulation depth and encapsulates all other parameters.

.. code-block:: python
    
    bg = lambda t,lam: bg_hom3d(t,200,lam)   # define function for background
    K = dipolarkernel(t,r,lam,bg)

The use of function handles is central to DeerLab, especially when fitting experimental data.


Time-domain signals
--------------------------------------------

To generate complete time-domain signals from a distance distribution and a background decay, use the function ``dipolarkernel`` and apply it to the distance distribution.

.. code-block:: python

        K = dipolarkernel(t,r,lam,B)   # generate dipolar kernel
        V = K@P                        # generate dipolar signal
        plt.plot(t,V)

It is possible to add noise to simulated data by using the ``whitegaussnoise`` function:

.. code-block:: python

        sigma = 0.05                        # noise level
        V = K@P + whitegaussnoise(t,sigma)  # add some noise

With this, uncorrelated Gaussian noise with standard deviation given as ``sigma`` is added to the noise-free signal.

Adding a phase rotation is also possible, yielding a complex-valued signal with non-zero imaginary component. The phase shift on the noise has to be taken into account too:

        phase = np.pi/4                   # phase shift, radians
        V = K@P*exp(-1j*phase)            # add a phase shift
        rnoise = whitegaussnoise(t,sigma) # real-component noise
        inoise = whitegaussnoise(t,sigma) # imaginary-component noise
        V = V + rnoise + inoise           # complex-valued noisy signal

