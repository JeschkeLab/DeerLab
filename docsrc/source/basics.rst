Basics
=========================================

DeerLab relies on a few central quantities. They include distance distributions and model functions, background and their model functions, time-domain signals, and kernels. They are described in this section.

All functions in DeerLab use the same units: all distances are in units of **nanometers**, and all times in units of **microseconds**.

Distance distributions
        A distance distribution `P(r)` between two spins is represented by a pair of vectors: a distance vector ``r`` (in nanometers) and a vector of densities ``P``
        (in inverse nanometers). The distance-vector ``r`` can be a linearly or non-linearly increasing vector, and it cannot have negative values. 
        The elements ``P[i]`` are the distance distribution values at ``r[i]`` and are all strictly non-negative. Outside of the range defined by ``r``, 
        the distribution ``P`` is assumed to be zero, i.e. the distribution is truncated to the range ``r``.
        The distance distribution ``P`` is normalized such that the integral over the range of the provided ``r`` equals one:

        .. math:: \int P(r) \mathrm{d}r = 1

        Furthermore, DeerLab distinguishes between *non-parametric* and *parametric* distance distributions.

        Non-parametric distance distributions
                These provide the more general definition of distance distributions. They have no particular shape and are represented by the vectors ``P`` and ``r``. 
                For example, you can generate ``P`` and ``r`` by an external program for spin label rotamer modeling. 
                
                In least-squares fitting, non-parametric distance distributions are preferred over parametric distance distributions, since
                they make fewer assumptions about the distribution and are more flexible. They introduce less bias.

        Parametric distance distributions 
                These have specific shapes that depend on a few parameters. DeerLab provides many parametric distance distribution :ref:`model functions <modelsref_dd>`. All 
                these functions start with the prefix ``dd_`` (``dd`` stands for "distance distribution"). They take a distance vector ``r`` and a parameter vector ``param``
                as inputs and return the distance distribution as a vector ``P``. Here is an example: ::

                        r = np.linspace(1.5,6,200) # distance-vector, in nm
                        rmean = 3 # nm, Gaussian mean distance
                        sigma = 0.2 # nm, Gaussian standard deviation
                        P = dl.dd_gauss(r,[rmean, sigma]) # Parametric Gaussian distribution
                        plt.plot(r,P)

                To programmatically get information on a particular distance distribution :ref:`model functions <modelsref_dd>` and its parameters, call the function without input arguments ::

                        info = dl.dd_gauss()        # obtain information on model and parameters

                which returns a dictionary ``info`` containing the the names of the model parameters, and other quantities such as their built-in lower and upper boundaries. 
                
                DeerLab provides a wide range of :ref:`parametric distribution models<modelsref_dd>` that fall into several groups. 


.. _bgmodels:

Dipolar background
        In DeerLab, all inter-molecular contributions to the dipolar signal (i.e. the signal due to randomly distributed spins in the sample that are not part of the 
        spin-labeled protein or object) are referred to as the dipolar background. There is a wide collection of parametric models that can be used to model the background. 
        All these background :ref:`model functions <modelsref_bg>` start with the prefix ``bg_``, they take the time axis vector ``t`` (in microseconds), and a parameter vector ``param``. 
        as inputs. The output is a background vector ``B`` defined over ``t``. To get information on the model and its parameters, call the function without inputs: ::

                info = dl.bg_hom3d()        # obtain information on model and parameters

        DeerLab's :ref:`background models<modelsref_bg>` fall into two categories, physical and phenomenological: 

        Physical background models
                Describe particular distributions of spin labels in space and depend on physical parameters such as spin concentration, exclusion distances, and fractal dimensionality.
                The most common background model is :ref:`bg_hom3d`, which describes the signal due to a homogeneous three-dimensional distribution of spins of a given concentration.
                A background due to a homogeneous distribution of spins in fractal dimensions is available with :ref:`bg_homfractal`, and excluded-volume effects can be accounted for using
                :ref:`bg_hom3dex` to model the background. 

                In addition to ``t`` and the model parameters ``param`` physical background model functions take the dipolar pathway amplitude ``lambda`` as a third input, for example ::

                        t = np.linspace(-0.1,4,200)    # time, in microseconds
                        lam = 0.4                      # modulation depth
                        conc = 70                      # spin concentration, in uM
                        B = dl.bg_hom3d(t,conc,lam)    # homogeneous 3D background
                        plt.plot(t,B)

        Phenomenological background models
                Represent various mathematical functions that are intended to *mimic* the background decay, without reference to a particular spatial distribution of spins. The parameters
                of these models do no have a direct physical meaning. Some examples include :ref:`bg_exp`, which models the background decay as a simple exponential function, or :ref:`bg_strexp`
                which model the background decay as a stretched exponential function.

                Phenomenological background model functions just take ``t`` and the model parameters ``param``  as input, for example ::

                        t = np.linspace(-0.1,4,200)    # time, in microseconds
                        kappa = 0.35                   # decay rate, in inverse microseconds
                        B = dl.bg_exp(t,kappa)         # exponential background
                        plt.plot(t,B) 
                
        In general, it is preferable to use the physical instead of phenomenological models.


.. _exmodels:

Experiments
        DeerLab supports a wide range of dipolar EPR experiments. Different experiments differ in their modulated dipolar pathways. Each of these pathways leads to a dipolar modulation 
        contribution to the dipolar signal, with a certain amplitude and refocusing times. A dipolar signal can be defined as a combination of an unmodulated contribution of and a contribution
        from all modulated pathways, which can de defined by their amplitude, refocusing time, and harmonic. For each type of supported dipolar EPR experiment, there is a dedicated experiment 
        :ref:` model function<modelsref_ex>` starting with ``ex_``, which models the dipolar pathways for that specific experiment. These functions take an array of parameters characterizing 
        the experiment. As output, they return an array containing information about the dipolar pathways of the experiment model.

        For example, the model function representing the typical model for a 4-pulse DEER signal is ``ex_4pdeer``: ::

                t = np.linspace(0,3,151)
                lam = 0.3;
                pathways = dl.ex_4pdeer(t,lam)

        The returned output ``pathways`` is a list of pathways information ::

                pathways = [[0.7], [0.3, 0]]

        Each nested list holds information about one pathway. The first element is modulation amplitude, and the second element is the refocusing point.
        In the above example, the first list shows a pathway with amplitude 0.7 and no refocusing time, indicating that it represents the unmodulated contribution.
        The pathway of the second list shows amplitude of 0.3 and refocusing time 0, i.e. this is the primary dipolar pathway.



Dipolar Kernel
        One of the core functions of DeerLab is ``dipolarkernel``. It constructs the kernel that provides the connection between the distance distribution and the time-domain dipolar signal 
        via the Fredholm integral equation 

        .. math:: V(t) = \int K(t,r)P(r) \mathrm{d}r

        The most simple dipolar kernel just requires the time-vector ``t`` and distance-vector ``r`` ::

                t = np.linspace(0,6,300)        # time axis, in us
                r = np.linspace(2,7,300)        # distance axis, in nm
                K0 = dl.dipolarkernel(t,r)      # dipolar kernel matrix

        To calculate the dipolar signal corresponding to a distance distribution ``P`` according to the equation above, use ::
        
                V = K0@P

        The above ``K0`` is the most elementary kernel, giving a dipolar signal without any background decay, and with a single dipolar evolution function centered at time zero and with modulation depth of 1.

        The kernel can also account for the background and the dipolar pathways. Then, operation  ``V=K@P`` will return the complete time-domain dipolar signal. Here is an example for a 4-pulse DEER signal ::

                lam = 0.4
                B = dl.bg_hom3d(t,200,lam)
                K = dl.dipolarkernel(t,r,lam,B)
                V = K@P
                plt.plot(t,V)

        When accounting for more than one dipolar pathway, the different refocusing times and modulation amplitudes must be provided to ``dipolarkernel``. Additionally, the background must be provided as a callable
        function that takes only time and modulation amplitude and encapsulates all other parameters. For example, for a 5-pulse DEER signal :: 

                Lam0 = 0.5      # amplitude of unmodulated component
                lam1 = 0.4      # amplitude of primary pathway
                lam2 = 0.1      # amplitude of secondary pathway
                T02 = 3.1       # refocusing time of secondary pathway, in us
                pathways = dl.ex_5pdeer([Lam0,lam1,lam2,T02]) # dipolar pathways of 5-pulse DEER experiment
                Bfcn = lambda t,lam: dl.bg_hom3d(t,200,lam)   # define function for background
                K = dl.dipolarkernel(t,r,pathways,Bfcn) # 5-pulse DEER dipolar kernel
        
        The function ``dipolarkernel`` also has :ref:`options<dipolarkernel>` to add an excitation bandwidth limitation, to select the internal calculation method, and more.


Dipolar signals
        Dipolar signals are the results of the many different dipolar EPR spectroscopy experiments. They represent the data from which distance distributions can be infered. 
        DeerLab provides the tools for simulating dipolar signals originating from different experiments. Note that these simulations are not based on spin dynamics simulations,
        but rather on theoretical analytical treatments of such problems.  

        Simulations in DeerLab can be easily achieved, as mentioned above, via a Fredholm integral using the correct dipolar kernel. To generate complete time-domain signals from 
        a distance distribution and a background decay, use the function ``dipolarkernel`` and apply it to the distance distribution. ::

                K = dl.dipolarkernel(t,r,lam,B)   # generate dipolar kernel
                V = K@P                           # generate dipolar signal
                plt.plot(t,V)

        It is possible to add noise to simulated data by using the ``whitegaussnoise`` function: ::

                sigma = 0.05                           # noise level
                V = K@P + dl.whitegaussnoise(t,sigma)  # add some noise

        With this, uncorrelated Gaussian noise with standard deviation given as ``sigma`` is added to the noise-free signal.

        Adding a phase rotation is also possible, yielding a complex-valued signal with non-zero imaginary component. The phase shift on the noise has to be taken into account too: ::

                phase = np.pi/4                      # phase shift, radians
                V = K@P*exp(-1j*phase)               # add a phase shift
                rnoise = dl.whitegaussnoise(t,sigma) # real-component noise
                inoise = dl.whitegaussnoise(t,sigma) # imaginary-component noise
                V = V + rnoise + inoise              # complex-valued noisy signal

