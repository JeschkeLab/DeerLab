import numpy as np
from deerlab.dd_models import dd_gauss
import inspect 
import matplotlib.pyplot as plt
import difflib



class FitResult(dict):
# ========================================================================
    r""" Represents the results of either the :ref:`fit` or :ref:`snlls` functions. 
    Depending on which function is used, the attributes and methods of the FitResult will be different.

 
    Attributes
    ----------
    model : ndarray
        The fitted model response.
    modelUncert : 
        Uncertainty quantification of the fitted model response.
    param : ndarray
        Fitted parameter vector ordered according to the model parameter indices.
    paramUncert : :ref:`UQResult`
        Uncertainty quantification of the parameter vector ordered according to the model parameter indices.

    regparam : float scalar
        Regularization parameter used in the fit.
    noiselvl: ndarray
        Estimated noise level of the data or user-provided noise level.
    success : bool
        Whether or not the optimizer exited successfully.
    cost : float
        Value of the cost function at the solution.
    residuals : ndarray
        Vector of residuals at the solution.
    stats : dict
        Goodness of fit statistical estimators:

        * ``stats['chi2red']`` - Reduced \chi^2 test
        * ``stats['r2']`` - R^2 test
        * ``stats['rmsd']`` - Root-mean squared deviation (RMSD)
        * ``stats['aic']`` - Akaike information criterion
        * ``stats['aicc']`` - Corrected Akaike information criterion
        * ``stats['bic']`` - Bayesian information criterion

    nonlin : ndarray
        Fitted non-linear parameters. [:ref:`snlls` specific attribute]
    nonlinUncert : :ref:`UQResult`
        Uncertainty quantification of the non-linear parameter set. [:ref:`snlls` specific attribute]
    lin : ndarray
        Fitted linear parameters. [:ref:`snlls` specific attribute]
    linUncert : :ref:`UQResult`
        Uncertainty quantification of the linear parameter set. [:ref:`snlls` specific attribute]

    """

    def __getattr__(self, attr):
        try:
            return self[attr]
        except KeyError:
            errstr = f"The results object has no attribute '{attr}'."
            attributes = [key for key in self.keys()]
            proposal = difflib.get_close_matches(attr, attributes)
            if len(proposal)>0:
                errstr += f' \n Did you mean: {proposal} ?'
            raise AttributeError(errstr)

    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__

    def __str__(self): 
        return self._summary

    def __repr__(self):
        if self.keys():
            m = max(map(len, list(self.keys()))) + 1
            return '\n'.join([k.rjust(m) + ': ' + repr(v)
                              for k, v in sorted(self.items())])
        else:
            return self.__class__.__name__ + "()"

    def __dir__(self):
        return list(self.keys())
    

    def _extarct_params_from_model(self, model):
        if callable(model):
            try:
                modelparam = model._parameter_list('vector')
                function_type=False
            except AttributeError:
                function_type=True
                modelparam = inspect.getfullargspec(model).args

        if not hasattr(self,'param'):
            raise ValueError('The fit object does not contain any fitted parameters.')

        # # Enforce model normalization
        # normfactor_keys = []
        # for key in modelparam:
        #     param = getattr(model,key)
        #     if np.all(param.linear):
        #         if param.normalization is not None:
        #             normfactor_key = f'{key}_scale'
        #             normfactor_keys.append(normfactor_key)
        #             try:
        #                 model.addnonlinear(normfactor_key,lb=-np.inf,ub=np.inf,par0=1,description=f'Normalization factor of {key}')
        #                 getattr(model,normfactor_key).freeze(1)
        #             except KeyError:
        #                 pass
                    

        # # Get some basic information on the parameter vector
        # modelparam = model._parameter_list(order='vector')
        # param_idx = [[] for _ in model._parameter_list('vector')]
        # idxprev = 0
        # for islinear in [False,True]:
        #     for n,param in enumerate(model._parameter_list('vector')):
        #         if np.all(getattr(model,param).linear == islinear):
        #             N = len(np.atleast_1d(getattr(model,param).idx))
        #             param_idx[n] = np.arange(idxprev,idxprev + N)
        #             idxprev += N  

        fitparams = {key : fitvalue if len(fitvalue)>1 else fitvalue[0] for key, fitvalue in zip(self.paramlist,[self.param[idx] for idx in self._param_idx])}
        # Check that all parameters are in the fit object
        for param in modelparam:
            if not param in fitparams: 
                raise KeyError(f'The fit object does not contain the {param} parameter.')
        

        params = {param : fitparams[param] for param in modelparam}
        params_idx = [self._param_idx[self.paramlist.index(param)] for param in modelparam]

        return modelparam, params, params_idx

    def _extract_params_from_function(self,function):
        """
        Extracts the fitted parameters from a callable function. 

        Assumes that all parameters are length 1.
        
        """
        # Get the parameter names from the function definition
        modelparam = inspect.getfullargspec(function).args
        
        fitparam_idx = self._param_idx
        
        # Get the parameter values from the fit object
        fitparams = {key : fitvalue if len(fitvalue)>1 else fitvalue[0] for key, fitvalue in zip(self.paramlist,[self.param[idx] for idx in fitparam_idx])}
        params = {param : fitparams[param] for param in modelparam}
        params_idx = [fitparam_idx[self.paramlist.index(param)] for param in modelparam]

        
    
            
        return modelparam, params, params_idx



    def evaluate(self, model, *constants):
    # ----------------------------------------------------------------------------
        """
        Evaluate a model at the fitted parameter values. 

        Takes a model object or callable function model to be evaluated. All the 
        parameters in the model or in the callable definition must match their 
        corresponding parameter names in the FitResult object. Any model 
        constants present required by the model must be specified as a second 
        argument constants. It returns the model's response at the fitted 
        parameter values as an ndarray.

        Parameters
        ----------

        model : :ref:`Model` or callable
            Model object or callable function to be evaluated. All the parameters in the model or in the callable definition
            must match their corresponding parameter names in the ``FitResult`` object.   
        constants : array_like 
            Any model constants present required by the model.  
        
        Returns
        -------

        response : array_like 
            Model response at the fitted parameter values. 
        """
        try:
            modelparam = model._parameter_list('vector')
            modelparam, fitparams, fitparam_idx = self._extarct_params_from_model(model)
        except AttributeError:
            modelparam, fitparams, fitparam_idx = self._extract_params_from_function(model)
        
        
        parameters = {param: fitparams[param] for param in modelparam}

        # Evaluate the input model
        response = model(*constants,**parameters)
        return response
    
    def propagate(self, model, *constants, lb=None, ub=None):
        """
        Propagate the uncertainty in the fit results to a model's response.

        Takes a model object or callable function model to be 
        evaluated. All the parameters in the model or in the callable definition 
        must match their corresponding parameter names in the FitResult object.
        Any model constants present required by the model must be specified as 
        a second argument constants. The lower bounds lb and upper bounds ub of
        the model's response can be specified as a third and fourth argument 
        respectively. It returns the model's response uncertainty 
        quantification as a UQResult object.
        
        Parameters
        ----------

        model : :ref:`Model` or callable
            Model object or callable function to be evaluated. All the parameters in the model or in the callable definition
            must match their corresponding parameter names in the ``FitResult`` object.   
        constants : array_like 
            Model constants. 
        lb : array_like, optional 
            Lower bounds of the model response.
        ub : array_like, optional 
            Upper bounds of the model response.   

        Returns
        -------

        responseUncert : :ref:`UQResult`
            Uncertainty quantification of the model's response.
        """

        try:
            modelparam = model._parameter_list('vector')
            modelparam, fitparams, fitparam_idx = self._extarct_params_from_model(model)

        except AttributeError:
            modelparam, fitparams, fitparam_idx = self._extract_params_from_function(model)


        # Propagate the uncertainty from that subset to the model
        modeluq = self.paramUncert.propagate(lambda param: model(*constants,*[param[s] for s in fitparam_idx]),lb,ub)
        return modeluq
    

    def plot(self,axis=None,xlabel=None,gof=False,fontsize=13):
            """
            Function to display the results. 
            
            This can also plot goodness-of-fit tests if requested.

            * Plot of residuals along with the estimated noise level and mean value
            * Histogram of the residuals weighted by the noise level, compared to the standard normal distribution
            * Autocorrelogram of the residuals, along the confidence region for a white noise vector
            
            Parameters
            ----------
            axis : array_like, optional
                Vector for the x-axis. If not specified, the default is the array indices.

            xlabel : str, optional
                Label for the x-axis. If not specified, the default is 'Array elements'.

            gof : bool, optional
                If set to True, the goodness-of-fit plots will be displayed, the default is False.
            
            fontsize : int, optional
                Fontsize of the figure labels, the default is 13.

            Returns
            -------
            fig : matplotlib.Figure
                Figure object containing the plot.


            """
                
            ys, yfits, yuqs, noiselvl, masks = getattr(self,'__plot_inputs')

            # Check which datasets are complex-valued
            complexy = [np.iscomplex(y).any() for y in ys]

            # Determine the distribution of the subplots in the figure
            nrows = len(ys) + np.sum(complexy)
            if gof:
                ncols = 4
                fig,axs = plt.subplots(nrows,ncols,figsize=[4*ncols,4*nrows], constrained_layout=True)
            else: 
                ncols = 1
                fig,axs = plt.subplots(nrows,ncols,figsize=[7*ncols,4*nrows])
            axs = np.atleast_1d(axs)
            axs = axs.flatten() 
            n = 0 # Index for subplots

            # If abscissa of the datasets are not specified, resort to default 
            if axis is None: 
                axis = [np.arange(len(y)) for y in ys]
            if not isinstance(axis,list): 
                axis = [axis]
            axis = [np.real(ax) for ax in axis]
            if xlabel is None: 
                xlabel = 'Array elements'

            # Go through every dataset
            for i,(y,yfit,yuq,noiselvl,mask) in enumerate(zip(ys,yfits,yuqs,noiselvl,masks)): 

                axis_masked = axis[i][mask]
                y_masked = y[mask]
                y_hidden = y[~mask]
                

                # If dataset is complex-valued, plot the real and imaginary parts separately
                if complexy[i]:
                    components = [np.real,np.imag]
                    componentstrs = [' (real)',' (imag)']
                else:
                    components = [np.real]
                    componentstrs = ['']
                for component,componentstr in zip(components,componentstrs):
                    
                    # Plot the experimental signal and fit
                    axs[n].plot(axis_masked,component(y_masked),'.',color='grey',label='Data'+componentstr)
                    if y_hidden.size>0:
                        axs[n].plot(axis[i][~mask],component(y_hidden),'.',color='grey',alpha=0.4, label='Masked data'+componentstr)
                    axs[n].plot(axis[i],component(yfit),color='#4550e6',label='Model fit')
                    if yuq.type!='void': 
                        axs[i].fill_between(axis[i],component(yuq.ci(95)[:,0]),component(yuq.ci(95)[:,1]),alpha=0.4,linewidth=0,color='#4550e6',label='95%-confidence interval')
                    axs[n].set_xlabel(xlabel,size=fontsize)
                    axs[n].set_ylabel(f'Dataset #{i+1}'+componentstr,size=fontsize)
                    axs[n].spines.right.set_visible(False)
                    axs[n].spines.top.set_visible(False) 
                    axs[n].legend(loc='best',frameon=False)
                    plt.autoscale(enable=True, axis='both', tight=True)
                    n += 1

                    # Plot the visual guides to assess the goodness-of-fit (if requested)
                    if gof: 
                        # Get the residual
                        residuals = component(yfit[mask] - y[mask])
                        

                        # Plot the residual values along the estimated noise level and mean value
                        axs[n].plot(axis_masked,residuals,'.',color='grey')
                        axs[n].hlines(np.mean(residuals),axis[i][0],axis_masked[-1],color='#4550e6',label='Mean')
                        axs[n].hlines(np.mean(residuals)+noiselvl,axis_masked[0],axis_masked[-1],color='#4550e6',linestyle='dashed',label='Estimated noise level')
                        axs[n].hlines(np.mean(residuals)-noiselvl,axis_masked[0],axis_masked[-1],color='#4550e6',linestyle='dashed')
                        axs[n].set_xlabel(xlabel,size=fontsize)        
                        axs[n].set_ylabel(f'Residual #{i+1}'+componentstr,size=fontsize)      
                        axs[n].spines.right.set_visible(False)
                        axs[n].spines.top.set_visible(False) 
                        axs[n].legend(loc='best',frameon=False)
                        plt.axis("tight")
                        n += 1

                        # Plot the histogram of the residuals weighted by the noise level, compared to the standard normal distribution
                        bins = np.linspace(-4,4,20)
                        axs[n].hist(residuals/noiselvl,bins,density=True,color='b',alpha=0.6, label='Residuals')
                        bins = np.linspace(-4,4,300)
                        N0 = dd_gauss(bins,0,1)
                        axs[n].get_yaxis().set_visible(False)
                        axs[n].fill(bins,N0,'k',alpha=0.4, label='$\mathcal{N}(0,1)$')
                        axs[n].set_xlabel('Normalized residuals',size=fontsize)       
                        axs[n].set_yticks([])
                        axs[n].spines.right.set_visible(False)
                        axs[n].spines.left.set_visible(False)
                        axs[n].spines.top.set_visible(False) 
                        axs[n].legend(loc='best',frameon=False)
                        n += 1

                        # Plot the autocorrelogram of the residuals, along the confidence region for a white noise vector
                        maxLag = len(residuals)-1
                        axs[n].acorr(residuals, usevlines=True, normed=True, maxlags=maxLag, lw=2,color='#4550e6',label='Residual autocorrelation')
                        threshold = 1.96/np.sqrt(len(residuals))
                        axs[n].fill_between(np.linspace(0,maxLag),-threshold,threshold,color='k',alpha=0.3,linewidth=0,label='White noise confidence region')
                        axs[n].get_yaxis().set_visible(False)
                        plt.axis("tight")
                        axs[n].set_xbound(lower=-0.5, upper=maxLag)
                        axs[n].spines.right.set_visible(False)
                        axs[n].spines.left.set_visible(False)
                        axs[n].spines.top.set_visible(False)
                        axs[n].set_xlabel('Lags',size=fontsize)       
                        axs[n].legend(loc='best',frameon=False)
                        n += 1

            # Adjust fontsize
            for ax in axs:
                for label in (ax.get_xticklabels() + ax.get_yticklabels()):
                    label.set_fontsize(fontsize)

            return fig
# ===========================================================================================
