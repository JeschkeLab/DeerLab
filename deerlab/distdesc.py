# distdesc.py - Distance distribution descriptors
# -----------------------------------------------
# This file is a part of DeerLab. License is MIT (see LICENSE.md).
# Copyright(c) 2019-2020: Luis Fabregas, Stefan Stoll and other contributors.

import numpy as np
import warnings
import copy
from scipy.signal import find_peaks

def distdesc(P,r,Puq=None,verbose=False):
    r""" Computes descriptors/estimators for the location, spread and shape 
    of a distance distribution with or without their corresponding uncertainties..

    Parameters
    ----------
    P : array_like
        Distance distribution.
    r : array_like
        Distance axis in nm. 
    Puq : :ref:`UncertQuant`
        Uncertainty quantification of the distance distribution. If not 
        specified, the one single output is returned without any uncertainty 
        estimation. If specified, two outputs are returned containing the 
        uncertainty estimation.

    Returns
    -------
    estimators : dict
        Dictionary of shape, location and spread descriptors of 
        the input distance distribution:
        
        Location parameters

            * ``'mean'`` - Mean distance in nm (see `more <https://en.wikipedia.org/wiki/Mean>`_)
            * ``'median'`` - Median distance in nm (see `more <https://en.wikipedia.org/wiki/Median>`_)
            * ``'iqm'`` - Interquartile mean (IQM) distance in nm (see `more <https://en.wikipedia.org/wiki/Interquartile_mean>`_)
            * ``'mode'`` - Modal distance in nm (see `more <https://en.wikipedia.org/wiki/Mode_(statistics)>`_)

        Spread parameters

            * ``'iqr'`` - Interquartile range in nm (see `more <https://en.wikipedia.org/wiki/Interquartile_range>`_)
            * ``'mad'`` - Mean absolute deviation (MAD)  in nm (see `more <https://en.wikipedia.org/wiki/Average_absolute_deviation>`_)
            * ``'std'`` - Standard deviation  in nm (see `more <https://en.wikipedia.org/wiki/Standard_deviation>`_)
            * ``'var'`` - Variance in nm² (see `more <https://en.wikipedia.org/wiki/Variance>`_)

        Shape parameters

            * ``'modality'`` - Modality (number of peaks)
            * ``'skewness'`` - Skewness (see `more <https://en.wikipedia.org/wiki/Skewness>`_)
            * ``'kurtosis'`` - Kurtosis (see `more <https://en.wikipedia.org/wiki/Kurtosis>`_)
    
    uq : dict of :ref:`UncertQuant`
        Dictionary of the parameters covariance-based uncertainty quantifications. 
        See above for the dictionary keys. Only returned if ``Puq`` is specified

    Other Parameters
    ----------------
    verbose : boolean
        Enable to print a summary of all distribution descriptors along 
        with their corresponding 95% confidence intervals. Disabled by default.

    Notes
    -----
    Both the ``'mode'`` and ``'modality'`` parameters have no corresponding covariance-based
    uncertainty quantification since they are mathematically not defined. These can, however, be 
    calculated via bootsrapping of these quantities, e.g. ::

        def analyze_rmode(V):
            fit = dl.fitsignal(V,t,r)
            rmode = dl.descdist(fit.P,r)['mode']
            return rmode
        # Bootstrap analysis of distance mode    
        rmode_uq = dl.bootan(V,Vfit,analyze_r_mode)

    """

    P,r = np.atleast_1d(P,r)

    # Auxiliary functions
    # -------------------

    # Percentile function
    def pctile(P,p):
        cdf = np.cumsum(P/np.sum(P))
        cdf, index = np.lib.arraysetops.unique(cdf,return_index=True)
        rpctile = np.interp(p/100,cdf,r[index])
        return rpctile
    # Expectation operator function
    def E(x,P):
        return np.sum(x*P/np.sum(P))
        
    # Location estimators
    # -------------------
    # 1st moment  - Mean 
    meanfcn = lambda P: E(r,P)
    # Median
    medianfcn = lambda P: pctile(P,50)
    # Interquartile mean
    iqmfcn = lambda P: E(r[(r>pctile(P,25)) & (r<pctile(P,75))],P[(r>pctile(P,25)) & (r<pctile(P,75))]) 
    # Mode
    modefcn = lambda P: r[np.argmax(P/np.sum(P))]

    # Spread estimators
    # -----------------
    # Interquartile range
    iqrfcn = lambda P: pctile(P,75) - pctile(P,25)
    # Mean absolute deviation
    madfcn = lambda P: E(abs(r - meanfcn(P)),P)
    # 2nd moment - Variance
    variancefcn = lambda P: E(r**2 - meanfcn(P)**2,P)
    # 2nd moment - Standard deviation
    stdfcn = lambda P: np.sqrt(variancefcn(P))

    # Shape estimators
    # ----------------
    # Modality
    modalityfcn = lambda P:  np.size(find_peaks(P, [np.max(P)/10])[0])
    # 3rd moment - Skewness
    skewnessfcn = lambda P: E(((r - meanfcn(P))/stdfcn(P))**3,P)
    # 4th moment - Kurtosis
    kurtosisfcn = lambda P: E(((r - meanfcn(P))/stdfcn(P))**4,P)

    # Ignore warnings when evaluating Jacobians outside of reasonable bounds
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')

        # Construct dict of distribution estimators
        estimators = {
            'mean': meanfcn(P),
            'median': medianfcn(P),
            'mode': modefcn(P),
            'iqm': iqmfcn(P),
            'mad': madfcn(P),
            'var': variancefcn(P),
            'std': stdfcn(P),
            'iqr': iqrfcn(P),
            'modality': modalityfcn(P),
            'skewness': skewnessfcn(P),
            'kurtosis': kurtosisfcn(P)
        }

        # If the uncertainty of the distance distribution is not provided....
        if Puq is None:            
            # Print summary of estimators if requested
            if verbose:
                print('-------------------------------------------------')
                print('Location')
                print('-------------------------------------------------')
                print('Range                    {:.2f}-{:.2f} nm'.format(min(r),max(r)))
                print('Mean                     {:.2f} nm'.format(estimators['mean']))
                print('Median                   {:.2f} nm'.format(estimators['median']))
                print('Interquartile mean       {:.2f} nm'.format(estimators['iqm']))
                print('Mode                     {:.2f} nm'.format(estimators['mode']))
                print('-------------------------------------------------')
                print('Spread')
                print('-------------------------------------------------')
                print('Standard deviation       {:.2f} nm'.format(estimators['std']))
                print('Mean absolute deviation  {:.2f} nm'.format(estimators['mad']))
                print('Interquartile range      {:.2f} nm'.format(estimators['iqr']))
                print('Variance                 {:.2f} nm²'.format(estimators['var']))
                print('-------------------------------------------------')
                print('Shape')
                print('-------------------------------------------------')
                print('Modality                 {0}    '.format(estimators['modality']))
                print('Skewness                 {:.2f} '.format(estimators['skewness']))
                print('Kurtosis                 {:.2f} '.format(estimators['kurtosis']))
                print('-------------------------------------------------')
            return estimators

        else:

            def propagation(Puq,fcn):
                # Propagate the uncertainty to the function
                uq_ = Puq.propagate(fcn)
                uq = copy.deepcopy(uq_)
                def ci(p):
                    paramci = uq_.ci(p)
                    return [paramci[:,0][0],paramci[:,1][0]]
                # Wrap the ci() method to simplify array
                uq.ci = ci
                return uq

            # Construct dict of their respective uncertainties
            uq = {
                'mean': propagation(Puq,meanfcn),
                'median': propagation(Puq,medianfcn),
                'mode': None,
                'iqm': propagation(Puq,iqmfcn),
                'mad': propagation(Puq,madfcn),
                'var': propagation(Puq,variancefcn),
                'std': propagation(Puq,stdfcn),
                'iqr': propagation(Puq,iqrfcn),
                'modality': None,
                'skewness': propagation(Puq,skewnessfcn),
                'kurtosis': propagation(Puq,kurtosisfcn)
            }

        # Print summary of estimators with 95%-confidence intervals if requested
        if verbose:
            print('-------------------------------------------------')
            print('Location')
            print('-------------------------------------------------')
            print('Range                    {:.2f}-{:.2f} nm'.format(min(r),max(r)))
            print('Mean                     {:.2f} ({:.2f},{:.2f}) nm'.format(estimators['mean'],uq['mean'].ci(95)[0],uq['mean'].ci(95)[1]))
            print('Median                   {:.2f} ({:.2f},{:.2f}) nm'.format(estimators['median'],uq['median'].ci(95)[0],uq['median'].ci(95)[1]))
            print('Interquartile mean       {:.2f} ({:.2f},{:.2f}) nm'.format(estimators['iqm'],uq['iqm'].ci(95)[0],uq['iqm'].ci(95)[1]))
            print('Mode                     {:.2f} nm'.format(estimators['mode']))
            print('-------------------------------------------------')
            print('Spread')
            print('-------------------------------------------------')
            print('Standard deviation       {:.2f} ({:.2f},{:.2f}) nm'.format(estimators['std'],uq['std'].ci(95)[0],uq['std'].ci(95)[1]))
            print('Mean absolute deviation  {:.2f} ({:.2f},{:.2f}) nm'.format(estimators['mad'],uq['mad'].ci(95)[0],uq['mad'].ci(95)[1]))
            print('Interquartile range      {:.2f} ({:.2f},{:.2f}) nm'.format(estimators['iqr'],uq['iqr'].ci(95)[0],uq['iqr'].ci(95)[1]))
            print('Variance                 {:.2f} ({:.2f},{:.2f}) nm²'.format(estimators['var'],uq['var'].ci(95)[0],uq['var'].ci(95)[1]))
            print('-------------------------------------------------')
            print('Shape')
            print('-------------------------------------------------')
            print('Modality                 {0}                    '.format(estimators['modality']))
            print('Skewness                 {:.2f} ({:.2f},{:.2f}) '.format(estimators['skewness'],uq['skewness'].ci(95)[0],uq['skewness'].ci(95)[1]))
            print('Kurtosis                 {:.2f} ({:.2f},{:.2f}) '.format(estimators['kurtosis'],uq['kurtosis'].ci(95)[0],uq['kurtosis'].ci(95)[1]))
            print('-------------------------------------------------')

        return estimators,uq
