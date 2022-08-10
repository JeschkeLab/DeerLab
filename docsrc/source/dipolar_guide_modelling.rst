.. _dipolar_modelling:

Modelling
=========================

DeerLab provides a very flexible framework to model dipolar signals originating from any dipolar EPR spectroscopy experiments. Choosing a model that properly describes your sample and experiment is of paramount importance. The DeerLab function ``dipolarmodel`` already defines the core model structure based on dipolar pathways, with the following components to be chosen:     

* **Distance range**: Also called the interspin distance axis, is the range of distances where the distribution is defined. 

* **Distribution model**: Describes the intra-molecular distance distribution in either a parametric (e.g. a Gaussian distribution) or a non-parametric way. 

* **Background model**: Describes the dipolar background signal arising from the inter-molecular contributions. 

* **Number of pathways**: Sets the number of dipolar pathways contributing to the dipolar signal.

For each of these four components, a choice needs to be made: 

Choosing a distance range
*************************

The distance range :math:`[r_\mathrm{min},r_\mathrm{max}]` is an important choice, as any distance distribution is truncated to this range, i.e. :math:`P(r)=0` for :math:`r<r_\mathrm{min}` and :math:`r>r_\mathrm{max}`. The lower limit of the distance range is determined by the bandwidth of the pulses, and also by the time increment. Typically, 1.5 nm is a reasonable choice. The upper limit depends on the distances in your sample. The number of points in ``r`` is usually set to a certain resolution (typically 0.01-0.05nm). Such a distance-axis is usually defined as ``r`` is most easily defined using the ``linspace`` function from NumPy: ::

    r = np.linspace(1.5,6.5,100)  # define distance range from 1.5nm to 6.5nm with a resolution of 0.05nm

Choosing a distribution model
******************************

A non-parametric distribution is specified by setting the choice of ``Pmodel`` keyword in ``dipolarmodel`` to ``None``. In a non-parametric distribution, each element :math:`P_i` of the distribution is a linear parameter. Non-parametric distributions are obtained via methods such as Tikhonov regularization. If there are reasons to believe that the distance distribution has a specific shape (e.g. Gaussian, Rice, random-coil, etc.), or if there is very little information in the data, use a parametric distance distribution model from the :ref:`list of available models<modelsref_dd>`.

Choose a background model
*************************

Typically, a background model of a homogenous 3D distribution of spins is appropriate. The associated parametric model function is :ref:`bg_hom3d`. In some cases, depending on the properties of your sample, other background models might be needed, such as backgrounds arising from distributions of spins in fractal dimensions or when  accounting for volume-exclusion effects. In such cases, use the associated parametric background models from the :ref:`list of available models<modelsref_bg>`. If there is no inter-molecular background in your sample, or it is negligible, set the background model to ``None``.



Choosing the number of dipolar pathways
*************************************** 

This decision should be based on the experiment you used to acquire the data and the type of pulses you used. In the case of 4-pulse DEER data, when analyzing a standard 4-pulse DEER signal without 2+1 component at the end a single pathway suffices. If the 2+1 component (appearing at the right edge of the time trace) is present, then it should be fitted as well, including its counterpart appearing at negative times, making a total of three dipolar pathways. Experiments such as 5-pulse DEER typically require at least two dipolar pathways to be properly modelled. 

Using experimental pulse delays
******************************** 

The dipolar pathways of a newly constructed dipolar model are initialized at arbitrary refocusing times and fully unconstrained. The refocusing times can be strongly constrained by knowing the experimental pulse sequence delays used to acquire the data. If the experiment used to acquire the data is known, as well as its pulse delays, then it is strongly recommended do so. 
 
DeerLab provides a selection of experimental information generators for some of the most widely employed experimental methods (see the of :ref:`list of available experiments <modelsref_ex>`). These are functions that take the pulse sequence delays, and return an ``ExperimentInfo`` object. This can be passed to the ``dipolarmodel`` function via the ``experiment`` keyword argument, to incorporate the experiment information on the model and constrain some of its parameters. 

When using experimental time delays and the ``experiment`` argument, the model assumes that the experimental time axis ``t`` has its zero time at the beginning of the interpulse delay (see the illustrations of the individual experiment models for details). However, experimentally it is common not to record the first few hundred nanoseconds of signal. This results in a so-called deadtime, which many commercial spectrometers (such as Bruker) do not account for when storing the time-vector. It is very important to account for that deadtime in the model via :: 

    deadtime = 0.4 # Experimental deadtime of 400ns, in μs
    t = t - deadtime # Shift the time axis to account for the deadtime 

If the time vector ``t`` does not have any deadtime, this step can be skipped. Otherwise, an incorrectly defined time vector will results in wrong results.

Constructing the dipolar model 
*******************************

Once all the decisions above have been made, the dipolar model can be constructed using the ``dipolarmodel`` function. The models that have an associated parametric function, e.g. ``bg_hom3d``, must be passed directly as inputs to ``dipolarmodel``. In Python, functions can be passed as inputs to other functions.  See the :ref:`details <dipolarmodel>` on ``dipolarmodel`` for more information. 

Example: Single-pathway 4-pulse DEER model
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
For example, a 4pDEER signal with non-parametric distance distribution and homogenous 3D background can be constructed using ::

    expinfo = dl.ex_4pdeer(tau1=0.5, tau2=5.5, pathways=[1])
    Vmodel = dl.dipolarmodel(t, r, Pmodel=None, Bmodel=dl.bg_hom3d, experiment=expinfo) 

By default, the function ``dipolarmodel`` assumes a non-parametric distance distribution, a homogenous 3D background and a single pathway. Thus the above is equivalent to ::

    expinfo = dl.ex_4pdeer(tau1=0.5, tau2=5.5, pathways=[1])
    Vmodel = dl.dipolarmodel(t, r, experiment=expinfo) 


Example: Two-pathway 5-pulse DEER model
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
For example, a 5pDEER signal with non-parametric distance distribution and homogenous 3D background can be constructed using ::

    expinfo = dl.ex_rev5pdeer(tau1=0.5, tau2=5.5, tau3=0.2, pathways=[1,2])
    Vmodel = dl.dipolarmodel(t, r, Pmodel=None, Bmodel=dl.bg_hom3d, experiment=expinfo)

Manipulating the model
***********************

A full summary of the constructed model(s) can be inspected by printing the model object ::

    >>> print(Vmodel)
    Description: Dipolar signal model
    Signature: (mod, reftime, conc, P)
    Constants: []
    Parameter Table: 
    ========= ======= ======= ======= ======== ======== ====== ====================================== 
    Name      Lower   Start   Upper    Type    Frozen   Unit   Description                           
    ========= ======= ======= ======= ======== ======== ====== ====================================== 
    mod           0    0.01       1   nonlin     No            Modulation depth                      
    reftime    -inf       0     inf   nonlin     No      μs    Refocusing time                       
    conc       0.01      50   5e+03   nonlin     No      μM    Spin concentration                    
    P             0       0     inf   linear     No     nm⁻¹   Non-parametric distance distribution  
    ========= ======= ======= ======= ======== ======== ====== ====================================== 


From this point on, the model can be modified, manipulated and expanded freely as any other DeerLab model. Check out the :ref:`modelling guide <modelling_guide>` for more details and instructions on model manipulation.

