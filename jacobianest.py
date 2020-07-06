import numpy as np
import scipy
import scipy.optimize as opt
import numpy.matlib

def jacobianest(fun,x0):
    """
     jacobianest: estimate of the Jacobian matrix of a vector valued function of n variables
     usage: [jac,err] = jacobianest(fun,x0)
    
     
     arguments: (input)
      fun - (vector valued) analytical function to differentiate.
            fun must be a function of the vector or array x0.
     
      x0  - vector location at which to differentiate fun
            If x0 is an nxm array, then fun is assumed to be
            a function of n*m variables.
    
    
     arguments: (output)
      jac - array of first partial derivatives of fun.
            Assuming that x0 is a vector of length p
            and fun returns a vector of length n, then
            jac will be an array of size (n,p)
    
      err - vector of error estimates corresponding to
            each partial derivative in jac.
    """

    # subfunction - swap vector elements
    # ----------------------------------------------------------------------
    def swapelement(vec,ind,val):
        """
        Swaps val as element ind, into the vector vec
        """
        vec[ind] = val
        return vec 
    # ----------------------------------------------------------------------

    # subfunction - romberg extrapolation
    # ----------------------------------------------------------------------
    def rombextrap(StepRatio,der_init,rombexpon):
        """
        Romberg extrapolation for each estimate
    
        StepRatio - Ratio decrease in step
        der_init - initial derivative estimates
        rombexpon - higher order terms to cancel using the romberg step
    
        der_romb - derivative estimates returned
        errest - error estimates
        amp - noise amplification factor due to the romberg step
        """
        srinv = 1/StepRatio
        # do nothing if no romberg terms
        nexpon = len(rombexpon)
        rmat = np.ones((nexpon+2,nexpon+1))
        # two romberg terms
        rmat[1,1:3] = np.power(srinv,rombexpon)
        rmat[2,1:3] = np.power(srinv,2*rombexpon)
        rmat[3,1:3] = np.power(srinv,3*rombexpon)
        # qr factorization used for the extrapolation as well
        # as the uncertainty estimates
        qromb,rromb = np.linalg.qr(rmat,'reduced')
        # the noise amplification is further amplified by the Romberg step.
        # amp = cond(rromb)
        # this does the extrapolation to a zero step size.
        ne = len(der_init)
        rhs = vec2mat(der_init,nexpon+2,ne - (nexpon+2))
        rombcoefs = scipy.linalg.solve(rromb,qromb.T@rhs)
        der_romb = rombcoefs[0,:].T
        # uncertainty estimate of derivative prediction
        s = np.sqrt(sum((rhs - rmat@rombcoefs)**2,0))
        rinv = scipy.linalg.solve(rromb,np.eye(nexpon+1))
        cov1 = np.sum(rinv**2,1) 
        errest = s.T*12.7062047361747*np.sqrt(cov1[0])

        return der_romb, errest 
    # ----------------------------------------------------------------------

    # subfunction - vec2mat
    # ----------------------------------------------------------------------
    def vec2mat(vec,n,m):
        # forms the matrix M, such that M(i,j) = vec(i+j-1)
        i,j = np.mgrid[range(n),range(-1,m)]
        ind = i+j
        mat = vec[ind]
        if n==1:
            mat = mat.T
        return mat
    # ----------------------------------------------------------------------
    
    nx = len(x0)
    MaxStep = 2
    StepRatio = 2.0000001

    # Get fun at the center point
    f0 = fun(x0)
    n = len(f0)
    if n==0:
        # Empty begets empty
        jac = np.zeros((0,nx))
        err = jac
        return jac, err
    relativedelta = MaxStep*StepRatio**np.arange(0,-25,-1)
    nsteps = len(relativedelta)

    # total number of derivatives we will need to take
    jac = np.zeros((n,nx))
    err = jac
    for i in range(nx):
        x0_i = x0[i]
        if x0_i != 0:
            delta = x0_i*relativedelta
        else:
            delta = relativedelta
        
        # evaluate at each step, centered around x0_i
        # difference to give a second order estimate
        fdel = np.zeros((n,nsteps))
        for j in range(nsteps):
            fdif = fun(swapelement(x0,i,x0_i + delta[j])) - fun(swapelement(x0,i,x0_i - delta[j]))
            fdel[:,j] = fdif[:]
        
        
        # these are pure second order estimates of the
        # first derivative, for each trial delta.
        derest = fdel*np.matlib.repmat(0.5/delta,n,1)
        
        # The error term on these estimates has a second order
        # component, but also some 4th and 6th order terms in it.
        # Use Romberg exrapolation to improve the estimates to
        # 6th order, as well as to provide the error estimate.
        
        # loop here, as rombextrap coupled with the trimming
        # will get complicated otherwise.
        for j in range(n):
            der_romb,errest = rombextrap(StepRatio,derest[j,:],np.asarray([2, 4]))
            
            # trim off 3 estimates at each end of the scale
            nest = len(der_romb)
            trim = [range(3), nest+(np.arange(-2,1))]
            der_romb = np.sort(der_romb)
            tags = np.argsort(der_romb)
            der_romb = np.delete(der_romb,trim)
            tags = np.delete(tags,trim)
            
            errest = errest[tags]
            
            # now pick the estimate with the lowest predicted error
            err[j,i] = np.min(errest)
            ind = np.argmin(errest)
            jac[j,i] = der_romb[ind]

    return jac, err 
# ----------------------------------------------------------------------
