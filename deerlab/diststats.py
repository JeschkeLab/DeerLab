# diststats.py - Distance distribution descriptors
# -----------------------------------------------
# This file is a part of DeerLab. License is MIT (see LICENSE.md).
# Copyright(c) 2019-2021: Luis Fabregas, Stefan Stoll and other contributors.

import numpy as np
import warnings
import copy
from scipy.signal import find_peaks
from scipy.integrate import cumtrapz

def diststats(r, P, Puq=None, verbose=False, threshold=None):
    r""" Computes statistical quantities for the location, spread, and shape 
    of a distance distribution, with or without their corresponding uncertainties.

    Parameters
    ----------
    r : array_like
        Distance axis, in nm. 
    
    P : array_like
        Distance distribution, does not have to be normalized.
    
    Puq : :ref:`UQResult`, optional
        Uncertainty quantification of the distance distribution. If Puq is not given, a
        single output is returned without any uncertainty estimation. If given, two outputs
        are returned containing the uncertainty estimation.
    
    verbose : boolean, optional
        Print a summary of all statistical quantities (and their uncertainties if calculated).
    
    threshold : float, optional
        Peak detection threshold for the calculation of modes of a distribution. The default is ``max(P)/10``.

    Returns
    -------
    estimators : dict
        Dictionary of shape, location, and spread descriptors of the input distance distribution:
        
        General parameters

            * ``'rmin'`` - Minimum distance in the distribution range, in nm
            * ``'rmax'`` - Maximum distance in the distribution range, in nm
            * ``'int'`` - Integral of the distance distribution

        Location parameters

            * ``'mean'`` or ``'moment1'`` - Mean distance, in nm (see `more <https://en.wikipedia.org/wiki/Mean>`__)
            * ``'median'`` - Median distance, in nm (see `more <https://en.wikipedia.org/wiki/Median>`__)
            * ``'iqm'`` - Interquartile mean (IQM) distance, in nm (see `more <https://en.wikipedia.org/wiki/Interquartile_mean>`__)
            * ``'mode'`` - First modal distance, in nm (see `more <https://en.wikipedia.org/wiki/Mode_(statistics)>`__)
            * ``'modes'`` - All modal distances, in nm (see `more <https://en.wikipedia.org/wiki/Mode_(statistics)>`__)

        Spread parameters

            * ``'iqr'`` - Interquartile range, in nm (see `more <https://en.wikipedia.org/wiki/Interquartile_range>`__)
            * ``'mad'`` - Mean absolute deviation (MAD), in nm (see `more <https://en.wikipedia.org/wiki/Average_absolute_deviation>`__)
            * ``'std'`` - Standard deviation, in nm (see `more <https://en.wikipedia.org/wiki/Standard_deviation>`__)
            * ``'var'`` or ``'moment2'`` - Variance, in nm² (see `more <https://en.wikipedia.org/wiki/Variance>`__)
            * ``'entropy'`` - Shannon entropy, in nat (see `more <https://en.wikipedia.org/wiki/Entropy_(information_theory)>`__)

        Shape parameters

            * ``'modality'`` - Modality (number of modes)
            * ``'skewness'`` or ``'moment3'`` - Skewness (see `more <https://en.wikipedia.org/wiki/Skewness>`__)
            * ``'kurtosis'`` - Excess kurtosis (see `more <https://en.wikipedia.org/wiki/Kurtosis>`__)
            * ``'moment4'`` - 4th moment (kurtosis) (see `more <https://en.wikipedia.org/wiki/Kurtosis>`__)

    uq : dict of :ref:`UQResult`
        Dictionary of the parameters covariance-based uncertainty quantifications. 
        See above for the dictionary keys. Only calculated if ``Puq`` is provided.

    Notes
    -----
    For the ``'mode'``, ``'modes'`` and ``'modality'`` parameters, covariance-based uncertainties are not 
    available. Unvertainties can, however, be calculated via bootsrapping of these quantities, e.g. ::

        def analyze_rmode(V):
            fit = dl.fitmodel(V,t,r)
            rmode = dl.diststats(fit.P,r)[0]['mode']
            return rmode
        # Bootstrap analysis of distance mode
        rmode_uq = dl.bootan(V,Vfit,analyze_rmode)

    """

    P,r = np.atleast_1d(P,r)

    # Check to avoid non-sensical results with syntax diststats(P,r)
    if not np.all(np.diff(r)>0): 
        raise KeyError('The distance axis (1st argument) must be a monotnously increasing vector.')

    if threshold is None:
        threshold = np.max(P)/10

    # Auxiliary functions
    # -------------------
    # Percentile function
    def pctile(r,P,p):
        cdf = cumtrapz(P,r,initial=0)
        cdf, index = np.lib.arraysetops.unique(cdf,return_index=True)
        rpctile = np.interp(p/100,cdf,r[index])
        return rpctile
    # Expectation operator function
    def E(x,P,r):
        return np.trapz(x*P,r)

    # Location estimators
    # -------------------
    # 1st moment  - Mean 
    meanfcn = lambda P: E(r,P,r)
    # Median
    medianfcn = lambda P: pctile(r,P,50)
    # Interquartile mean
    def iqmfcn(P):
        IQrange = (r>pctile(r,P,25)) & (r<pctile(r,P,75))
        return E(r[IQrange],P[IQrange]/np.trapz(P[IQrange],r[IQrange]),r[IQrange]) 
    # Mode
    modefcn = lambda P: r[np.argmax(P)]
    # Modes
    modesfcn = lambda P: r[find_peaks(P,height=threshold)[0]]

    # Spread estimators
    # -----------------
    # Interquartile range
    iqrfcn = lambda P: pctile(r,P,75) - pctile(r,P,25)
    # Mean absolute deviation
    madfcn = lambda P: E(abs(r - meanfcn(P)),P,r)
    # 2nd moment - Variance
    variancefcn = lambda P: E(r**2 - meanfcn(P)**2,P,r)
    # 2nd moment - Standard deviation
    stdfcn = lambda P: np.sqrt(variancefcn(P))
    # Entropy (information theory)
    entropyfcn = lambda P: -E(np.log(np.maximum(np.finfo(float).eps,P)),P,r)

    # Shape estimators
    # ----------------
    # Modality
    modalityfcn = lambda P:  np.size(modesfcn(P))
    # 3rd moment - Skewness
    skewnessfcn = lambda P: E(((r - meanfcn(P))/stdfcn(P))**3,P,r)
    # 4th moment - Kurtosis
    kurtosisfcn = lambda P: E(((r - meanfcn(P))/stdfcn(P))**4,P,r)
    # Excess kurtosis 
    exkurtosisfcn = lambda P: 3 - E(((r - meanfcn(P))/stdfcn(P))**4,P,r)
        
    # 0th moment  - Integral 
    intfcn = lambda P: np.trapz(P,r)

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
        print(f'Range                    {min(r):.2f}-{max(r):.2f} nm')
        print(f'Integral                 {estimators["int"]:.2f}')
        print('-------------------------------------------------')
        print('Location')
        print('-------------------------------------------------')
        print(f'Mean                     {estimators["mean"]:.2f} nm')
        print(f'Median                   {estimators["median"]:.2f} nm')
        print(f'Interquartile mean       {estimators["iqm"]:.2f} nm')
        print(f'Mode                     {estimators["mode"]:.2f} nm')
        print('-------------------------------------------------')
        print('Spread')
        print('-------------------------------------------------')
        print(f'Standard deviation       {estimators["std"]:.2f} nm')
        print(f'Mean absolute deviation  {estimators["mad"]:.2f} nm')
        print(f'Interquartile range      {estimators["iqr"]:.2f} nm')
        print(f'Variance                 {estimators["var"]:.2f} nm²')
        print('-------------------------------------------------')
        print('Shape')
        print('-------------------------------------------------')
        print(f'Modality                 {estimators["modality"]}    ')
        print(f'Skewness                 {estimators["skewness"]:.2f} ')
        print(f'Excess kurtosis          {estimators["kurtosis"]:.2f} ')
        print('-------------------------------------------------')
    else:
        print('-------------------------------------------------')
        print('Distribution Statistics')
        print('-------------------------------------------------')
        print(f'Range                    {min(r):.2f}-{max(r):.2f} nm')
        print(f'Integral                 {estimators["int"]:.2f}')
        print('-------------------------------------------------')
        print('Location')
        print('-------------------------------------------------')
        print(f'Range                    {min(r):.2f}-{max(r):.2f} nm')
        print(f'Mean                     {estimators["mean"]:.2f} ({uq["mean"].ci(95)[0]:.2f},{uq["mean"].ci(95)[1]:.2f}) nm')
        print(f'Median                   {estimators["median"]:.2f} ({uq["median"].ci(95)[0]:.2f},{uq["median"].ci(95)[1]:.2f}) nm')
        print(f'Interquartile mean       {estimators["iqm"]:.2f} ({uq["iqm"].ci(95)[0]:.2f},{uq["iqm"].ci(95)[1]:.2f}) nm')
        print(f'Mode                     {estimators["mode"]:.2f} nm')
        print('-------------------------------------------------')
        print('Spread')
        print('-------------------------------------------------')
        print(f'Standard deviation       {estimators["std"]:.2f} ({uq["std"].ci(95)[0]:.2f},{uq["std"].ci(95)[1]:.2f}) nm')
        print(f'Mean absolute deviation  {estimators["mad"]:.2f} ({uq["mad"].ci(95)[0]:.2f},{uq["mad"].ci(95)[1]:.2f}) nm')
        print(f'Interquartile range      {estimators["iqr"]:.2f} ({uq["iqr"].ci(95)[0]:.2f},{uq["iqr"].ci(95)[1]:.2f}) nm')
        print(f'Variance                 {estimators["var"]:.2f} ({uq["var"].ci(95)[0]:.2f},{uq["var"].ci(95)[1]:.2f}) nm²')
        print('-------------------------------------------------')
        print('Shape')
        print('-------------------------------------------------')
        print(f'Modality                 {estimators["modality"]}')
        print(f'Skewness                 {estimators["skewness"]:.2f} ({uq["skewness"].ci(95)[0]:.2f},{uq["skewness"].ci(95)[1]:.2f}) ')
        print(f'Kurtosis                 {estimators["kurtosis"]:.2f} ({uq["kurtosis"].ci(95)[0]:.2f},{uq["kurtosis"].ci(95)[1]:.2f}) ')
        print('-------------------------------------------------')
