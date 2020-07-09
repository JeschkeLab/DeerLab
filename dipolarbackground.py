import numpy as np
import math as m
import types
import pandas as pd

def dipolarbackground(
    t,pathinfo,Bmodel,
    renormalize = True,
    overtonecoeff = 1
    ):

    """
    DIPOLARBACKGROUND Multipathway background generator
    
       B = DIPOLARBACKGROUND(t)
       B = DIPOLARBACKGROUND(t,pathinfo,Bmodel)
       B = DIPOLARBACKGROUND(___,'Property',value,___)
    
    Computes the N-point multipathway background B from a N-point time axis t
    (in microseconds), the dipolar pathways defined in pathinfo and a basis 
    function Bmodel.
    
      Inputs:
         t         N-element time axis, in microseconds
         pathinfo  px2 or px3 array of modulation depths lambda, refocusing points
                   T0, and harmonics n for multiple pathways, each row contains
                   [lambda T0 n] or [lambda T0] for one pathway. If n is not given
                   it is assumed to be 1.
         Bmodel    Function handle for a background model: @(t)bg_model(t,par)
    
      Returns:
         B         N-element total multipathway background 

    This file is a part of DeerLab. License is MIT (see LICENSE.md).
    Copyright(c) 2019-2020: Luis Fabregas, Stefan Stoll and other contributors.
    """

    # Ensure that all inputs are numpy arrays
    t = np.atleast_1d(t)
    pathinfo = np.atleast_1d(pathinfo)
    overtonecoeff = np.atleast_1d(overtonecoeff)


    if type(Bmodel) is not types.LambdaType:
        raise TypeError('For a model with multiple modulated pathways, B must be a function handle of the type: @(t,lambda) bg_model(t,par,lambda)')

    if not np.isreal(pathinfo).all:
        raise TypeError('lambda/pathinfo must be a numeric array.')

    if len(pathinfo)==1:
        lam = pathinfo
        pathinfo = np.array([[1-lam, np.NaN], [lam, 0]], dtype=object)

    if not any(np.shape(pathinfo)[1]!=np.array([2, 3])):
        raise TypeError('pathinfo must be a numeric array with two or three columns.')

    if any(pd.isnull(pathinfo[:,0])):
        raise ValueError('In pathinfo, NaN can only appear in the second column (refocusing time) e.g. path[1,:] = [Lam0 NaN]')

    # Normalize the pathway amplitudes to unity
    pathinfo[:,0] = pathinfo[:,0]/sum(pathinfo[:,0])
    lam = pathinfo[:,0]
    T0 = pathinfo[:,1]
    if np.shape(pathinfo)[1]==2:
        n = np.ones(np.shape(T0))
    else:
        n = pathinfo[:,2]
    

    # Combine all unmodulated components, and eliminate from list
    unmodulated = np.where(pd.isnull(T0))
    lam = np.delete(lam,unmodulated)
    T0 = np.delete(T0,unmodulated)
    n = np.delete(n,unmodulated)  

    # Fold overtones into pathway list
    nCoeffs = len(overtonecoeff)
    lam_,T0_,n_ = (np.empty(0) for _ in range(3))
    for i in range(nCoeffs):
        lam_ = np.concatenate((lam_, lam*overtonecoeff[i]))
        T0_ = np.concatenate((T0_,T0))
        n_ = np.concatenate((n_,n*(i+1)))
    lam,T0,n = (lam_,T0_,n_)
    nModPathways = len(lam)

    # Construction of multi-pathway background function 
    #-------------------------------------------------------------------------------
    Bnorm = 1
    B = 1
    for pathway in range(nModPathways):        
        B = B*Bmodel(n[pathway]*(t-T0[pathway]),lam[pathway])
        Bnorm = Bnorm*Bmodel(-T0[pathway]*n[pathway],lam[pathway])
    
    if renormalize:
        B = B/Bnorm

    return B
