.. only:: html

    .. note::
        :class: sphx-glr-download-link-note

        Click :ref:`here <sphx_glr_download_auto_examples_plot_pseudotitration_parameter_free.py>`     to download the full example code
    .. rst-class:: sphx-glr-example-title

    .. _sphx_glr_auto_examples_plot_pseudotitration_parameter_free.py:


Analyzing pseudo-titration (dose-respononse) curves with parameter-free distributions 
======================================================================================

How to use separable non-linear least squares (SNLLS)
to fit a pseudo-titration curve to multiple DEER datsets, using
parameter-free distance distributions.


.. code-block:: python


    import numpy as np
    import matplotlib.pyplot as plt
    from deerlab import *








Generating multiple datasets
-----------------------------------------------------------------------------
 First, let's prepare the chemical desciption of the problem. In this example we will
 simulate a protein system in their states A (natural) and B (changed upon addition
 of a ligand L) given by the chemical equilibrium  A + L <-> B.


.. code-block:: python


    def chemicalequilibrium(Kdis,L):
        """Prepare equilibrium of type: A + L <-> B"""
        Ctot = 1 # total protein concentration, uM

        # # Get fraction of state B
        Kb = 1/Kdis
        xB = np.zeros_like(L)
        for q in range(len(L)):
            xB_ = np.roots([Kb*Ctot, -(Kb*L[q] + Kb*Ctot + 1), Kb*L[q]])
            try:
                xB[q] = xB_[(xB_<=1) & (xB_>=0)]
            except:
                xB[q] = np.minimum(1,np.maximum(0,xB_[0]))    
        # Get fraction of state A
        xA = 1 - xB

        return xA,xB








Next, we define the dipolar kernel model as the non-linear function of the
SNLLS problem. This function needs to take the parameters and return a
cell-array of kernels, each one for the corresponding datasets that we
have. 
Since we have a total distribution of the form 
    ``P = xA*PA + xB*PB``
we can define an augmented kernel as
    ``K = [xA*KA xB*KB]``
such that 
    ``K@[PA PB] = V``
and the vector ``[PA PB]`` constitutes the linear part fitted by SNLLS.


.. code-block:: python


    def Kmodel(par,ts,rA,rB,L):

        Nsignals = len(ts)

        # Unpack parameters
        lam,k,Kdis = par

        # Get fractions for given KD
        [xA,xB] = chemicalequilibrium(Kdis,L)

        Ks = [[]]*Nsignals
        # General the dipolar kernels
        for i in range(Nsignals):
            B = bg_exp(ts[i],k,lam)
            # Kernel for fraction A
            KstateA = dipolarkernel(ts[i],rA,lam,B)
            # Kernel for fraction B
            KstateB = dipolarkernel(ts[i],rB,lam,B)
            Ks[i] = np.concatenate((xA[i]*KstateA, xB[i]*KstateB),axis=1)

        return Ks








Now, we can simulate multiple signals corresponding to different concentrations
of added ligand. 


.. code-block:: python


    # Time axes
    ts = [[]]*7
    ts[0] = np.linspace(-0.2,3,100)
    ts[1] = np.linspace(-0.1,5,300)
    ts[2] = np.linspace(-0.5,2,200)
    ts[3] = np.linspace(-0.1,1,100)
    ts[4] = np.linspace(-0.2,6,300)
    ts[5] = np.linspace(-0.2,3,300)
    ts[6] = np.linspace(-0.1,4,100)
    Nsignals = len(ts)

    # Distance axes for states A and B
    rA = np.linspace(1,8,100)
    rB = np.linspace(1,8,100)

    # Distributions for states A and B
    PstateA = dd_gauss(rA,[5.5, 0.4])
    PstateB = dd_gauss2(rB,[4.5, 0.7, 0.4, 3.5, 0.6, 0.6])

    L = [0.3, 1, 3, 10, 30, 100, 300] # total ligand concentration, uM
    Kdis = 5.65  # dissociation constant, uM

    # Populations of states A and B
    [xA,xB] = chemicalequilibrium(Kdis,L)

    # Global kernel model
    Ks = Kmodel([0.25, 0.1, Kdis],ts,rA,rB,L)

    # Simulate dipolar signals
    Vs = [[]]*Nsignals
    for i in range(Nsignals):
        np.random.seed(i)
        Vs[i] = Ks[i]@np.concatenate((PstateA, PstateB)) + whitegaussnoise(ts[i],0.005)








Psuedotitration SNLLS Analysis
-----------------------------------------------------------------------------
 For simplification, we will assume that all DEER traces have the same
 background function and modulation depth. Thus, we will fit the
 modulations depth (lam) and background decay constant (k) globally along
 the dissociation constant (KD).


.. code-block:: python


    # Non-linear parameters:
    #       lam  k   KD
    par0 = [0.5, 0.5,  5]  # start values 
    lb   = [ 0,   0,   1]  # lower bounds
    ub   = [ 1,   1,  10] # upper bounds

    # Linear parameters:
    #     |-------PA--------||--------PB--------|
    lbl = np.concatenate((np.zeros_like(rA), np.zeros_like(rB))) # Non-negativity constraint
    ubl = [] # Unconstrained

    # Run SNLLS optimization
    fit = snlls(Vs,lambda p: Kmodel(p,ts,rA,rB,L),par0,lb,ub,lbl,ubl)
    # Extract fit results
    parfit = fit.nonlin
    Pfit = fit.lin
    puq = fit.uncertainty

    # Extract the fitted disociation constant value and its 95#-confidence interval
    Kdisfit = parfit[2]
    parci = puq.ci(95,'nonlin')
    KDci = parci[2,:]

    # Print result
    print('Kdis = {:.2f}({:.2f}-{:.2f})uM'.format(Kdisfit,KDci[0],KDci[1]))





.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    D:\lufa\projects\DeerLab\DeerLab\examples\plot_pseudotitration_parameter_free.py:34: ComplexWarning: Casting complex values to real discards the imaginary part
      xB[q] = np.minimum(1,np.maximum(0,xB_[0]))
    Kdis = 5.33(5.04-5.63)uM




Plot results


.. code-block:: python

    plt.figure(figsize=(12,12))

    # Simulate fits
    Ksfit = Kmodel(parfit,ts,rA,rB,L)
    Vsfit = []
    plt.subplot(3,2,(1,3))
    for i in range(Nsignals):
        Vsfit.append(Ksfit[i]@Pfit)
        plt.plot(ts[i],Vs[i]+i/9,'k.',ts[i],Vsfit[i]+i/9,'r',linewidth=1.5)
    plt.grid(alpha =0.3)
    plt.xlabel('t [$\mu s$]')
    plt.ylabel('V(t) [a.u.]')
    plt.legend(['data','fit'])

    xAfit,xBfit = chemicalequilibrium(Kdisfit,L)
    plt.subplot(2,2,(2,4))
    for i in range(Nsignals):
        PAfit = xAfit[i]*Pfit[0:len(rA)]
        PBfit = xBfit[i]*Pfit[len(rA):len(rB)+len(rA)]
        plt.plot(rA,PAfit+2*i,'r',rB,PBfit+2*i,'b',linewidth=1.5)

    plt.grid(alpha =0.3)
    plt.xlabel('r [nm]')
    plt.ylabel('P(r)')
    plt.legend(['state A','state B'])
    plt.xlim([2,7])

    plt.subplot(325)
    plt.plot(np.log10(L),xA,'r-',np.log10(L),xB,'b-')
    plt.plot(np.log10(L),xAfit,'ro',np.log10(L),xBfit,'bo',linewidth=1.5)
    plt.grid(alpha =0.3)
    plt.xlabel('log$_{10}$([L])')
    plt.ylabel('Fractions')
    plt.legend(['state A','state B'])
    plt.ylim([0,1])






.. image:: /auto_examples/images/sphx_glr_plot_pseudotitration_parameter_free_001.png
    :alt: plot pseudotitration parameter free
    :class: sphx-glr-single-img


.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none


    (0.0, 1.0)




.. rst-class:: sphx-glr-timing

   **Total running time of the script:** ( 0 minutes  26.262 seconds)


.. _sphx_glr_download_auto_examples_plot_pseudotitration_parameter_free.py:


.. only :: html

 .. container:: sphx-glr-footer
    :class: sphx-glr-footer-example



  .. container:: sphx-glr-download sphx-glr-download-python

     :download:`Download Python source code: plot_pseudotitration_parameter_free.py <plot_pseudotitration_parameter_free.py>`



  .. container:: sphx-glr-download sphx-glr-download-jupyter

     :download:`Download Jupyter notebook: plot_pseudotitration_parameter_free.ipynb <plot_pseudotitration_parameter_free.ipynb>`


.. only:: html

 .. rst-class:: sphx-glr-signature

    `Gallery generated by Sphinx-Gallery <https://sphinx-gallery.github.io>`_
