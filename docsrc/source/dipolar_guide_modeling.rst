.. _dipolar_modeling:

Modeling
=========================

DeerLab provides a very flexible framework to simulate/model signals from a wide range dipolar EPR spectroscopy experiments. Use the DeerLab function ``dipolarmodel`` to define the model.

    Vmodel = dl.dipolarmodel(t,r,Pmodel,Bmodel,experiment)

The following components need to be chosen:

* **Time range** ``t``: This is the range of dipolar evolution times for the desired signal.

* **Distance range** ``r``: This is the range of intra-molecular spin-spin distances where the distribution is defined.

* **Distribution model** ``Pmodel``: Describes the intra-molecular distance distribution in either a parametric form (e.g. a Gaussian distribution) or a non-parametric form (as a histogram). 

* **Background model** ``Bmodel``: Describes the dipolar background signal arising from the inter-molecular spin pairs. 

* **Type of experiment** ``experiment``: This sets the number of dipolar pathways contributing to the dipolar signal.

For each of these components, a choice needs to be made: 

Choosing a time range
*************************

The time range depends on the expected range of distances, but is typically in the range of 0.5 to several microseconds. To construct a time vector, use the ``linspace`` function from NumPy: ::

    t = np.linspace(0.2,3.0,321)  # time from 02. µm to 3.0 µm with a resolution of 0.01 µs

Note that DeerLab interprets ``t`` such that its zero point is right after the end of the preceding pulse (e.g. right after the second observer pulse in 4-pulse DEER) rather than at the refocusing point of the signal (after `\tau_1` after the second observer pulse).

Choosing a distance range
*************************

The distance range is a vector that goes from a minimum distance :math:`r_\mathrm{min}` to a maximum distance :math:`r_\mathrm{max}`. Any distance distribution is truncated to this range, i.e. :math:`P(r)=0` for :math:`r<r_\mathrm{min}` and :math:`r>r_\mathrm{max}`. The lower limit of the distance range is determined by the bandwidth of the pulses, and also by the time increment. Typically, 1.5 nm is a reasonable choice. The upper limit depends on the distances in your sample. The number of points in ``r`` is usually set to a certain resolution (typically 0.01-0.05 nm). To construct a distance vector, use the ``linspace`` function from NumPy: ::

    r = np.linspace(1.5,6.5,100)  # define distance range from 1.5 nm to 6.5 nm with a resolution of 0.05 nm

Choosing a distribution model
******************************

A non-parametric distribution is specified by setting the choice of ``Pmodel`` keyword in ``dipolarmodel`` to ``None``. In a non-parametric distribution, `P` is a vector of the same length as ``r``, and each element :math:`P_i` is a parameter. Basically, the distribution is represented as a histogram. When analyzing data, such distributions are used with methods such as Tikhonov regularization. If there are reasons to believe that the distance distribution has a specific shape (e.g. Gaussian, Rice, random-coil, etc.), or if there is very little information in the data, use a parametric distance distribution model from the :ref:`list of available models<modelsref_dd>`.

Choose a background model
*************************

Typically, a background model of a homogenous 3D distribution of spins is appropriate. The associated parametric model function is :ref:`bg_hom3d`. In some cases, depending on the properties of your sample, other background models might be needed, such as backgrounds arising from distributions of spins in fractal dimensions or when  accounting for volume-exclusion effects. In such cases, use the associated parametric background models from the :ref:`list of available models<modelsref_bg>`. If there is no inter-molecular background in your sample, or it is negligible, set the background model to ``None``.


Choosing the number of dipolar pathways
*************************************** 

This decision should be based on the experiment you used to acquire the data and the type of pulses you used. In the case of 4-pulse DEER data, when analyzing a standard 4-pulse DEER signal without 2+1 component at the end a single pathway suffices. If the 2+1 component (appearing at the right edge of the time trace) is present, then it should be fitted as well, including its counterpart appearing at negative times, making a total of three dipolar pathways. Experiments such as 5-pulse DEER typically require at least two dipolar pathways to be properly modelled.

Using experimental pulse delays
******************************** 

The dipolar pathways of a newly constructed dipolar model are initialized at arbitrary refocusing times and fully unconstrained. The refocusing times can be strongly constrained by knowing the experimental pulse sequence delays used to acquire the data. If the experiment used to acquire the data is known, as well as its pulse delays, then it is strongly recommended do so. 
 
DeerLab provides a selection of experimental information generators for some of the most widely employed experimental methods (see the of :ref:`list of available experiments <modelsref_ex>`). These are functions that take the pulse sequence delays, and return an ``ExperimentInfo`` object. This can be passed to the ``dipolarmodel`` function via the ``experiment`` keyword argument, to incorporate the experiment information on the model and constrain some of its parameters. These experimental information generators can also take information on the duration of the longest microwave pulses to more accurately constraint the parameters when using long pulses such as in frequency-swept or shaped microwave pulses.

When using an experimental time axis ``t`` and the ``experiment`` argument, the model assumes that ``t`` is zero at the beginning of the interpulse delay (see the illustrations of the individual experiment models for details). For example, for 4-pulse DEER, this zero point would be immediately after the second pulse. However, experimentally, it is common that the time axis is defined and stored differently. For 4-pulse DEER, a commercial spectrometer usually defines ``t`` such that it is zero ``tau1`` after the second pulse. It is important to correct for this offset (also called start time or deadtime):: 

    t0 = 0.4   # Experimental time offset of 400 ns, in μs
    # To convert from experimental time axis to DeerLab time axis
    t = t + t0 # Shift the experimental time axis to match DeerLab's definition
    
    # To convert from DeerLab time axis to experimental time axis
    t = t - t0
    
Without this shift, an incorrectly defined time vector will result in wrong simulated or fitted signals, and in wrong plots.


Constructing the dipolar model 
*******************************

Once all the decisions above have been made, the dipolar model can be constructed using the ``dipolarmodel`` function. The models that have an associated parametric function, e.g. ``bg_hom3d``, must be passed directly as inputs to ``dipolarmodel``. In Python, functions can be passed as inputs to other functions.  See the :ref:`details <dipolarmodel>` on ``dipolarmodel`` for more information. 

For example, a model for 4-pulse DEER with a non-parametric distance distribution and a homogenous 3D background can be constructed using ::

    expinfo = dl.ex_4pdeer(tau1=0.5, tau2=5.5, pathways=[1])
    Vmodel = dl.dipolarmodel(t, r, Pmodel=None, Bmodel=dl.bg_hom3d, experiment=expinfo) 

By default, the function ``dipolarmodel`` assumes a non-parametric distance distribution, a homogenous 3D background and a single pathway. Thus the above is equivalent to ::

    expinfo = dl.ex_4pdeer(tau1=0.5, tau2=5.5, pathways=[1])
    Vmodel = dl.dipolarmodel(t, r, experiment=expinfo) 


To construct a model for 5-pulse DEER with non-parametric distance distribution and homogenous 3D background can be constructed using ::

    expinfo = dl.ex_rev5pdeer(tau1=0.5, tau2=5.5, tau3=0.2, pathways=[1,2])
    Vmodel = dl.dipolarmodel(t, r, Pmodel=None, Bmodel=dl.bg_hom3d, experiment=expinfo)

Manipulating the model
***********************

To obtain a summary of the constructed model(s), print it: ::

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


Once defined, the model can be modified, manipulated and expanded freely as any other DeerLab model. Refer to the :ref:`modeling guide <modeling_guide>` for more details and instructions on model manipulation.
