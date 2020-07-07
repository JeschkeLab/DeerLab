import math as m
import numpy as np



def multigaussfun(r,r0,fwhm,a):
    
    n = len(r0)
    P = np.zeros_like(r)
    for k in range(0,n):
        sig = fwhm[k]/2/m.sqrt(2*m.log(2))
        P += a[k]*m.sqrt(1/(2*m.pi))*1/sig*np.exp(-0.5*((r-r0[k])/sig)**2)
    if not all(P==0):
        P = P/np.trapz(P,r);    
    
    return P

def dd_gauss(*args):

    if not args:
        info = dict(
            Parameters = ('Center','FWHM'),
            Units = ('nm','nm'),
            Start = np.asarray([3.4, 0.5]),
            Lower = np.asarray([1, 0.2]),
            Upper = np.asarray([20, 5])
        )
        return info
    elif len(args)==2:
        r = args[0]
        p = args[1]
    else:
        raise TypeError('Two input arguments required dd_model(r,p)')

    r0 = [p[0]]
    fwhm = [p[1]]
    a = [1.0]
    P = multigaussfun(r,r0,fwhm,a)
    return P
    
def dd_gauss2(*args):
    
    if not args:
        info = dict(
            Parameters = ('Center of 1st Gaussian', 'FWHM of 1st Gaussian', 'Amplitude of 1st Gaussian',
                          'Center of 2nd Gaussian', 'FWHM of 2nd Gaussian', 'Amplitude of 2nd Gaussian'),
            Units = ('nm','nm','','nm','nm',''),
            Start = np.asarray([2.5, 0.5, 0.5, 3.5, 0.5, 0.5]),
            Lower = np.asarray([1, 0.2, 0, 1, 0.2, 0]),
            Upper = np.asarray([20, 5, 1, 20, 5, 1])
        )
        return info
    elif len(args)==2:
        r = args[0]
        p = args[1]
    else:
        raise TypeError('Two input arguments required dd_model(r,p)')

    r0 = [p[0], p[3]]
    fwhm = [p[1], p[4]]
    a = [p[2], p[5]]
    P = multigaussfun(r,r0,fwhm,a)
    return P
    
def dd_gauss3(*args):
    
    if not args:
        info = dict(
            Parameters = ('Center of 1st Gaussian', 'FWHM of 1st Gaussian', 'Amplitude of 1st Gaussian',
                          'Center of 2nd Gaussian', 'FWHM of 2nd Gaussian', 'Amplitude of 2nd Gaussian',
                          'Center of 3rd Gaussian', 'FWHM of 3rd Gaussian', 'Amplitude of 3rd Gaussian'),
            Units = ('nm','nm','','nm','nm','','nm','nm',''),
            Start = np.asarray([2.5, 0.5, 0.3, 3.5, 0.5, 0.3, 5, 0.5, 0.3]),
            Lower = np.asarray([1, 0.2, 0, 1, 0.2, 0, 1, 0.2, 0]),
            Upper = np.asarray([20, 5, 1, 20, 5, 1,  20, 5, 1])
        )
        return info
    elif len(args)==2:
        r = args[0]
        p = args[1]
    else:
        raise TypeError('Two input arguments required dd_model(r,p)')

    r0 = [p[0], p[3], p[6]]
    fwhm = [p[1], p[4], p[7]]
    a = [p[2], p[5], p[8]]
    P = multigaussfun(r,r0,fwhm,a)
    return P

