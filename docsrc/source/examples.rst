Example scripts
=========================================

Here is a collection of example scripts for the use of DeerLab. You can also find them in the ``examples`` subfolder in the DeerLab folder.

.. note:: Couldn't find what you were looking for? `Add a request <https://github.com/JeschkeLab/DeerLab/issues/112>`_ for a new example/tutorial and it will considered for the next release.

.. toctree:
    :maxdepth: 0
    :hidden:
    :glob:
    :caption: Examples - Basics

    ./examples/example_tikhonovbasic
    ./examples/example_multigauss_4pdeer
    ./examples/example_mixmodels_fitting
    ./examples/example_ci_regularization
    ./examples/example_pakepattern
    ./examples/example_aptsimple
    ./examples/example_visualization_lcurve
    ./examples/example_backgroundtreatment
    ./examples/example_basicprocessing
    ./examples/example_timedomainfitting
    ./examples/example_selectmodel
    ./examples/example_globalfit_regularization
    ./examples/example_globalfit_localglobal_vars


.. raw:: html

    <div >
      <ul class="examples">
        <li onclick="window.location='examples/example_basicprocessing.html#example_basicprocessing'">
          <img src="_images/example_basicprocessing.png" >
          <h3>Fitting a 4-pulse DEER signal with a parameter-free distribution</h3>
          <p> How to fit a simple 4-pulse DEER signal with a parameter-free distribution, a background, and a modulation amplitude.</p>
        </li>

        <li onclick="window.location='examples/example_basic_5pdeer.html#example_basic_5pdeer'">
          <img src="_images/example_basic_5pdeer.png" >
          <h3>Fitting a 5-pulse DEER signal with a parameter-free distribution</h3>
          <p> How to fit a 5-pulse DEER signal with a parameter-free distribution, a background, and all multi-pathway parameters.</p>
        </li>
             
        <li onclick="window.location='examples/example_multigauss_4pdeer.html#example_multigauss_4pdeer'">
          <img src="_images/example_multigauss_4pdeer.png" >
          <h3>Multi-Gauss fitting of a 4-pulse DEER signal</h3>
          <p> How to fit a simple 4-pulse DEER signal with background using a multi-Gauss model, i.e automatically optimizing the number of Gaussians in the model.</p>
        </li>
        
        <li onclick="window.location='examples/example_mixmodels_fitting.html#example_mixmodels_fitting'">
          <img src="_images/example_mixmodels_fitting.png" >
          <h3>Fitting a mixed distance-domain model</h3>
          <p> Basic manipulation of parametric models and creating mixed models for fitting distance distributions.</p>
        </li>
        
        <li onclick="window.location='examples/example_emulating_deeranalysis.html#example_emulating_deeranalysis'">
          <img src="_images/example_emulating_deeranalysis.png" >
          <h3>Emulating the DeerAnalysis workflow</h3>
          <p> How to reproduce the type of workflow implemented in DeerAnalysis, using DeerLab functions. This kind of analysis workflow is outdated and not recommended for routine or accurate data analysis.</p>
        </li>
        
        <li onclick="window.location='examples/example_ci_regularization.html#example_ci_regularization'">
          <img src="_images/example_ci_regularization.png" >
          <h3>Comparing confidence intervals for regularization results</h3>
          <p> A simpe example of uncertainty estimation for Tikhonov regularization results. The example covers the use of confidence intervals obtained from curvature matrices and boostrap analysis. </p>
        </li>
    
        <li onclick="window.location='examples/example_bootan_pardist.html#example_bootan_pardist'">
          <img src="_images/example_bootan_pardist.png" >
          <h3>Bootstrapped distributions of fit parameters</h3>
          <p>  How to generate probability density functions of values for fit parameters using bootstrapping, showcased for a 5-pulse DEER signal.</p>
        </li>
    
        <li onclick="window.location='examples/example_param_error_propagation.html#example_param_error_propagation'">
          <img src="_images/example_param_error_propagation.png" >
          <h3>Uncertainty propagation from parameter fits</h3>
          <p> How to propagate the uncertainty of the fitted parameters to the models which depend on them. </p>
        </li>

        <li onclick="window.location='examples/example_pakepattern.html#example_pakepattern'">
          <img src="_images/example_pakepattern.png" >
          <h3>Analyzing the Pake pattern of a dipolar signal</h3>
          <p> A very basic example for displaying the frequency-domain Pake pattern (spectrum) of a given dipolar signal. </p>
        </li>

        <li onclick="window.location='examples/example_aptsimple.html#example_aptsimple'">
          <img src="_images/example_aptsimple.png" >
          <h3>Data analysis with the approximate Pake transformation (APT)</h3>
          <p> How to use the old approximate Pake transformation (APT) technique, which was commonly found in the old DeerAnalysis. </p>
        </li>

        <li onclick="window.location='examples/example_backgroundtreatment.html#example_backgroundtreatment'">
          <img src="_images/example_backgroundtreatment.png" >
          <h3>Comparison of background treatment approaches</h3>
          <p> A comparison of different approaches for treating the background in dipolar signals</p>
        </li>
        
        <li onclick="window.location='examples/example_visualization_lcurve.html#example_visualization_lcurve'">
          <img src="_images/example_visualization_lcurve.png" >
          <h3>Visualization of the regularization parameter selection on the L-curve</h3>
          <p> How to construct the L-curve for visualization of the optimal regularization parameter selection in a similar fashion to the old DeerAnalysis. </p>
        </li>
        
                <li onclick="window.location='examples/example_extracting_gauss_constraints.html#example_extracting_gauss_constraints'">
          <img src="_images/example_extracting_gauss_constraints.png" >
          <h3>Extracting Gaussian constraints from a parameter-free distribution fit</h3>
          <p> How to extract Gaussian constraints from a parameter-free fit of a dipolar signal and how to estimate the corresponding uncertainty. </p>
        </li>
                
        <li onclick="window.location='examples/example_customkernel_snlls.html#example_customkernel_snlls'">
          <img src="_images/example_customkernel_snlls.png" >
          <h3>Fitting a custom dipolar kernel model with a parameter-free distribution</h3>
          <p>How to use of SNLLS to fit a custom dipolar kernel model and a parameter-free distribution to a dipolar signal</p>
        </li>
                
        <li onclick="window.location='examples/example_parfree_titration.html#example_parfree_titration'">
          <img src="_images/example_parfree_titration.png" >
          <h3>Analyzing pseudo-titration (dose-respononse) curves with parameter-free distributions</h3>
          <p> How to use separable non-linear least squares (SNLLS) to fit a pseudo-titration curve to multiple DEER datsets, using parameter-free distance distributions.</p>
        </li>
                
        <li onclick="window.location='examples/example_timedomainfitting.html#example_timedomainfitting'">
          <img src="_images/example_timedomainfitting.png" >
          <h3>Fitting a custom time-domain model of a 4-pulse DEER signal</h3>
          <p> Hot to construct a custom time-domain parametric model and fit it using fitparamodel </p>
        </li>
     
        <li onclick="window.location='examples/example_selectmodel.html#example_selectmodel'">
          <img src="_images/example_selectmodel.png" >
          <h3>Selecting an optimal parametric model for fitting a dipolar signal</h3>
          <p> How to optimally select a parametric model for a given dipolar signal from a given set.  </p>
        </li>
                
        <li onclick="window.location='examples/example_globalfit_regularization.html#example_globalfit_regularization'">
          <img src="_images/example_globalfit_regularization.png" >
          <h3>Global fit of dipolar evolution functions using fitregmodel</h3>
          <p> How to do global fitting using Tikhonov regularization via fitregmodel. </p>
        </li>

        <li onclick="window.location='examples/example_globalfit_localglobal_vars.html#example_globalfit_localglobal_vars'">
          <img src="_images/example_globalfit_localglobal_vars.png" >
          <h3>Global model fits with global, local and fixed parameters</h3>
          <p>How to fit multiple signals to a global model, which may depend on some parameters which need to be globally fitted, some locally and some might be fixed and not fitted.  </p>
        </li>      
        
        
      </ul>
    </div>

