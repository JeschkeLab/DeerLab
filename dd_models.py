import math as m
import numpy as np
import matrixtools as mat

def dd_gauss(r,p):
    r0 = [p[0]]
    fwhm = [p[1]]
    a = [1.0]
    P = multigaussfun(r,r0,fwhm,a)
    P = mat.column(P)
    return P
    
def dd_gauss2(r,p):
    r0 = [p[0], p[3]]
    fwhm = [p[1], p[4]]
    a = [p[2], p[5]]
    P = multigaussfun(r,r0,fwhm,a)
    P = mat.column(P)
    return P
    
def dd_gauss3(r,p):
    r0 = [p[0], p[3], p[6]]
    fwhm = [p[1], p[4], p[7]]
    a = [p[2], p[5], p[8]]
    P = multigaussfun(r,r0,fwhm,a)
    P = mat.column(P)
    return P

def multigaussfun(r,r0,fwhm,a):
    n = len(r0)
    P = np.zeros_like(r)
    for k in range(0,n):
        sig = fwhm[k]/2/m.sqrt(2*m.log(2))
        P += a[k]*m.sqrt(1/(2*m.pi))*1/sig*np.exp(-0.5*((r-r0[k])/sig)**2)
    return P
