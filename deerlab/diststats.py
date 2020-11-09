# diststats.py - Distance distribution descriptors
# -----------------------------------------------
# This file is a part of DeerLab. License is MIT (see LICENSE.md).
# Copyright(c) 2019-2020: Luis Fabregas, Stefan Stoll and other contributors.

import numpy as np
import warnings
import copy
from scipy.signal import find_peaks

def diststats(r, P, Puq=None, verbose=False, threshold=None):
    r""" Computes statistical quantities for the location, spread and shape 
    of a distance distribution with or without their corresponding uncertainties.

    Parameters
    ----------
    r : array_like
        Distance axis in nm. 
    
    P : array_like
        Distance distribution.
    
    Puq : :ref:`UncertQuant`
        Uncertainty quantification of the distance distribution. If not 
        specified, the one single output is returned without any uncertainty 
        estimation. If specified, two outputs are returned containing the 
        uncertainty estimation.
    
    verbose : boolean, optional
        Enables printing a summary of all statistical quantities and their uncertainties if calculated.
    
    threshold : float, optional
        Peak detection threshold for the calculation of modes of a distribution. The default is ``max(P)/10``.

    Returns
    -------
    estimators : dict
        Dictionary of shape, location and spread descriptors of 
        the input distance distribution:
        
        General parameters

            * ``'rmin'`` - Minimum distance in the distribution range in nm
            * ``'rmax'`` - Maximum distance in the distribution range in nm
            * ``'int'`` - Integral of the distance distribution

        Location parameters

            * ``'mean'`` or ``'moment1'`` - Mean distance in nm (see `more <https://en.wikipedia.org/wiki/Mean>`_)
            * ``'median'`` - Median distance in nm (see `more <https://en.wikipedia.org/wiki/Median>`_)
            * ``'iqm'`` - Interquartile mean (IQM) distance in nm (see `more <https://en.wikipedia.org/wiki/Interquartile_mean>`_)
            * ``'mode'`` - First modal distance in nm (see `more <https://en.wikipedia.org/wiki/Mode_(statistics)>`_)
            * ``'modes'`` - All modal distances in nm (see `more <https://en.wikipedia.org/wiki/Mode_(statistics)>`_)

        Spread parameters

            * ``'iqr'`` - Interquartile range in nm (see `more <https://en.wikipedia.org/wiki/Interquartile_range>`_)
            * ``'mad'`` - Mean absolute deviation (MAD)  in nm (see `more <https://en.wikipedia.org/wiki/Average_absolute_deviation>`_)
            * ``'std'`` - Standard deviation  in nm (see `more <https://en.wikipedia.org/wiki/Standard_deviation>`_)
            * ``'var'`` or ``'moment2'`` - Variance in nm² (see `more <https://en.wikipedia.org/wiki/Variance>`_)
            * ``'entropy'`` - Shannon entropy (see `more <https://en.wikipedia.org/wiki/Entropy_(information_theory)>`_)

        Shape parameters

            * ``'modality'`` - Modality (number of peaks)
            * ``'skewness'`` or ``'moment3'`` - Skewness (see `more <https://en.wikipedia.org/wiki/Skewness>`_)
            * ``'kurtosis'`` - Excess kurtosis (see `more <https://en.wikipedia.org/wiki/Kurtosis>`_)
            * ``'moment4'`` - 4th moment (kurtosis) (see `more <https://en.wikipedia.org/wiki/Kurtosis>`_)

    uq : dict of :ref:`UncertQuant`
        Dictionary of the parameters covariance-based uncertainty quantifications. 
        See above for the dictionary keys. Only calculated if ``Puq`` is specified.

    Notes
    -----
    The ``'mode'``, ``'modes'`` and ``'modality'`` parameters have no corresponding covariance-based
    uncertainty quantification since they are mathematically not defined. These can, however, be 
    calculated via bootsrapping of these quantities, e.g. ::

        def analyze_rmode(V):
            fit = dl.fitsignal(V,t,r)
            rmode = dl.diststats(fit.P,r)[0]['mode']
            return rmode
        # Bootstrap analysis of distance mode    
        rmode_uq = dl.bootan(V,Vfit,analyze_r_mode)

    """

    P,r = np.atleast_1d(P,r)

    if threshold is None:
        threshold = np.max(P)/10

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
        
    # 0th moment  - Integral 
    intfcn = lambda P: np.trapz(P,r)

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
    # Modes
    modesfcn = lambda P: r[find_peaks(P,height=threshold)[0]]

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
    # Entropy (information theory)
    entropyfcn = lambda P: -E(np.log(np.maximum(np.finfo(float).eps,P)),P)

    # Shape estimators
    # ----------------
    # Modality
    modalityfcn = lambda P:  np.size(modesfcn(P))
    # 3rd moment - Skewness
    skewnessfcn = lambda P: E(((r - meanfcn(P))/stdfcn(P))**3,P)
    # 4th moment - Kurtosis
    kurtosisfcn = lambda P: E(((r - meanfcn(P))/stdfcn(P))**4,P)
    # Excess kurtosis 
    exkurtosisfcn = lambda P: 3 - E(((r - meanfcn(P))/stdfcn(P))**4,P)

    # Calculate distribution estimators
    estimators = {
        'rmin': min(r),
        'rmax': max(r),
        'int': intfcn(P),
        'mean': meanfcn(P),
        'median': medianfcn(P),
        'mode': modefcn(P),
        'modes': modesfcn(P),
        'iqm': iqmfcn(P),
        'mad': madfcn(P),
        'var': variancefcn(P),
        'std': stdfcn(P),
        'iqr': iqrfcn(P),
        'entropy': entropyfcn(P),
        'modality': modalityfcn(P),
        'skewness': skewnessfcn(P),
        'kurtosis': exkurtosisfcn(P),
        'moment1': meanfcn(P),
        'moment2': variancefcn(P),
        'moment3': skewnessfcn(P),
        'moment4': kurtosisfcn(P),
    }

    # Calculate distribution estimator uncertainties if uncertainty of P is given
    uq = None
    if Puq is not None:
        uq = {
            'rmin': None,
            'rmax': None,
            'int': None,
            'mean': _propagation(Puq,meanfcn),
            'median': _propagation(Puq,medianfcn),
            'mode': None,
            'modes': None,
            'iqm': _propagation(Puq,iqmfcn),
            'mad': _propagation(Puq,madfcn),
            'var': _propagation(Puq,variancefcn),
            'std': _propagation(Puq,stdfcn),
            'iqr': _propagation(Puq,iqrfcn),
            'entropy': _propagation(Puq,entropyfcn),
            'modality': None,
            'skewness': _propagation(Puq,skewnessfcn),
            'kurtosis': _propagation(Puq,exkurtosisfcn),
            'moment1': _propagation(Puq,meanfcn),
            'moment2': _propagation(Puq,variancefcn),
            'moment3': _propagation(Puq,skewnessfcn),
            'moment4': _propagation(Puq,kurtosisfcn),
        }

    # Print if requested
    if verbose:
        _print_estimators(r,estimators,uq)
        
    return estimators,uq


def _propagation(Puq,fcn):
    # Propagate the uncertainty to the function
    with warnings.catch_warnings():
        # Ignore warnings when evaluating Jacobians outside of reasonable bounds
        warnings.simplefilter('ignore')
        uq_ = Puq.propagate(fcn)
    
    uq = copy.deepcopy(uq_)
    def ci(p):
        paramci = uq_.ci(p)
        return [paramci[:,0][0],paramci[:,1][0]]
    # Wrap the ci() method to simplify array
    uq.ci = ci
    return uq


def _print_estimators(r,estimators,uq):
    # Print summary of estimators if requested
    if uq is None:
        
        print('-------------------------------------------------')
        print('Distribution Statistics')
        print('-------------------------------------------------')
        print('Range                    {:.2f}-{:.2f} nm'.format(min(r),max(r)))
        print('Integral                 {:.2f}'.format(estimators['int']))
        print('-------------------------------------------------')
        print('Location')
        print('-------------------------------------------------')
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
        print('Excess kurtosis          {:.2f} '.format(estimators['kurtosis']))
        print('-------------------------------------------------')
    else:
        print('-------------------------------------------------')
        print('Distribution Statistics')
        print('-------------------------------------------------')
        print('Range                    {:.2f}-{:.2f} nm'.format(min(r),max(r)))
        print('Integral                 {:.2f}'.format(estimators['int']))
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
