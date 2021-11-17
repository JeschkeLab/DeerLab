.. _dipolar_guide_fitting.rst

Fitting
============


Fitting the model
*****************

The model ``Vmodel`` can be fitted to the experimental data ``V`` by calling the ``fit`` function: ::

    result = dl.fit(Vmodel,Vexp)  # Fit the model to the experimental data


After ``fit`` has found a solution, it returns an object that we assigned to ``result``. This object contains fields with all quantities of interest with the fit results, such as the fitted model and parameters, goodness-of-fit statistics, and uncertainty information. Check out the :ref:`fitting guide <fitting_fitresult>` for more details on the quantities provided in ``result``.


Adding penalties
*****************

Penalty terms can be added to the objective function to impose certain properties upon the solution. While DeerLab can take any kind of penalty function (see the :ref:`fitting guide <fitting_guide>` for details), for dipolar models it provides a specialized function ``dipolarpenalty`` which easily generates penalties based on the distance distribution. 

To generate such a penalty, you must provide the model ``Pmodel`` for the distance distribution (as provided in ``dipolarmodel``), as well as the distance axis vector ``r``. Next, the type of penalty must be specified: 

- ``'compactness'``: Imposes compactness of the distance distribution. A compact distribution avoid having distribution mass spread towards the edges of the distance axis vector. 
- ``'smoothness'``: Imposes smoothness of the distance distribution. This is particularly useful for imposing smoothness of parametric models of the distance distribution. For non-parametric distributions, smoothness is already imposed by the regularization criterion, making this penalty unnecessary. 

All penalties are weighted by a weighting parameter, which is optimized according to a selection criterion which must be specified to the ``dipolarpenalty method``. For the ``smoothness`` penalty, the ``'aic'`` criterion is recommended, while for the ``smoothness`` criterion, the ``'icc'`` criterion is recommended.

The ``dipolarpenalty`` function will return a ``Penalty`` object which can be passed to the fit function through the ``penalties`` keyword argument. 


Example: Fitting a non-parametric distribution with a compactness criterion
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For example, to introduce compactness in the fit of a dipolar model with a non-parametric distance distribution we must set the distribution model to ``None`` to indicate a non-parametric distribution ::

    compactness_penalty = dl.dipolarpenalty(None, r, 'compactness', 'icc')
    results = dl.fit(Vmodel,Vexp, penalties=compactness_penalty)

Example: Fitting a Gaussian distribution with a compactness criterion
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For example, to introduce compactness in the fit of a dipolar model with a Gaussian distance distribution we must set the distribution model to ``dd_gauss`` to indicate the parametric distribution ::

    compactness_penalty = dl.dipolarpenalty(dl.dd_gauss, r, 'compactness', 'icc')
    results = dl.fit(Vmodel,Vexp, penalties=compactness_penalty)

Displaying the results
**********************

For just a quick display of the results, you can use the ``plot()`` method of the ``fit`` object that will display a figure with you experimental data, the corresponding fit including confidence bands. :: 

    results.plot(axis=t,xlabel='Time (μs)') # display results


.. image:: ./images/beginners_guide1.png
   :width: 450px

For a quick summary of the fit results, including goodness-of-fit statistics and the fitted model parameter values (including 95% confidence intervals), can be accessed by just printing the ``results`` object :: 

    >>>print(results)
    Goodness-of-fit: 
    ========= ============= ============ ========== ========== 
    Dataset   Noise level   Reduced 𝛘2     RMSD       AIC     
    ========= ============= ============ ========== ========== 
       #1       1426.905       1.036      1443.706   3631.484  
    ========= ============= ============ ========== ========== 
    Model parameters: 
    =========== ========= ========================= ======= ====================================== 
     Parameter   Value     95%-Confidence interval   Units   Description                           
    =========== ========= ========================= ======= ====================================== 
     mod         0.505     (0.494,0.516)                     Modulation depth                      
     reftime     0.096     (0.092,0.100)              μs     Refocusing time                       
     conc        295.909   (279.412,312.405)          μM     Spin concentration                    
     P           ...       (...,...)                 None    Non-parametric distance distribution  
    =========== ========= ========================= ======= ====================================== 

Any specific quantities can be extracted from the ``results`` object. For each parameter in the model, the ``results`` output contains an attribute ``results.<parameter>`` named after the parameter containing the fitted value of that parameter, as well as another attribute ``results.<parameter>Uncert`` containing the uncertainty estimates of that parameter, from which confidence intervals can be constructed (the :ref:`uncertainty guide <uncertainty>` for details). For example: :: 

    # Distance distribution 
    results.P # Fitted distance distribution 
    results.PUncert.ci(95) # Distance distribution 95% confidence intervals

    # Modulation depth 
    results.mod # Fitted modulation depth 
    results.modUncert.ci(95) # Modulation depth 95% confidence intervals



Exporting the figure and the data
*********************************

After completing the fit, you might want to export the figure with the fit. Here is one way to do it: ::

    figure = fit.plot()                       # get figure object
    figure.savefig('DEERFig.png', dpi=600)    # save figure as png file
    figure.savefig('DEERFig.pdf')             # save figure as pdf file

To export the fitted distance distribution for plotting with another software, save it in a simple text file ::

    np.savetxt('distancedistribution.txt', np.asarray((r, fit.P, *fit.Puncert.ci(95).T)).T)

The generated file contain four columns: the distance axis, the distance distributions, and the upper and lower confidence bounds. The ``.T`` indicate array transposes, which are used to get the confidence bands into the column format for saving.

To export the fitted time-domain trace, use similarly ::

    np.savetxt('timetrace.txt', np.asarray((t, V, fit.V, *fit.Vuncert.ci(95).T)).T)

