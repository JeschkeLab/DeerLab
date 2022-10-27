
import numpy as np
from numpy import pi, inf, NaN
from deerlab.bg_models import bg_hom3d,bg_exp
from deerlab.dd_models import dd_gauss
from deerlab.dipolarkernel import dipolarkernel,elementarykernel_twospin
from deerlab.utils import assert_docstring
from deerlab.constants import ge

def test_matrixsize_fresnel():
#=======================================================================
    "Check non-square kernel matrix construction"
    Ntime = 50
    Ndist = 70
    t = np.linspace(-0.5,5,50)
    r = np.linspace(2,6,70)
    K = dipolarkernel(t,r)
    assert K.shape == (Ntime, Ndist)
#=======================================================================


def test_singletime():
#=======================================================================
    "Check that one can generate a kernel with one single time-domain point"

    t = 0.5
    r = np.linspace(1,5,100)
    K = dipolarkernel(t,r)

    assert np.shape(K)[0]==1
#=======================================================================
    

def test_singledist():
#=======================================================================
    "Check that one can generate a kernel with one single distance-domain point"

    t = np.linspace(0,5,50)
    r = 3.5
    K = dipolarkernel(t,r)

    assert np.shape(K)[1]==1
#=======================================================================


def test_complex_type_fresnel():
#=======================================================================
    "Test that the values are complex-valued if requested."

    # Generate kernel numerically
    t = 1 # µs
    r = 1 # nm
    K = dipolarkernel(t,r,method='fresnel',complex=True)

    assert isinstance(K[0,0],np.complex128)
#=======================================================================

def test_complex_type_grid():
#=======================================================================
    "Test that the values are complex-valued if requested."

    # Generate kernel numerically
    t = 1 # µs
    r = 1 # nm
    K = dipolarkernel(t,r,method='grid',complex=True)

    assert isinstance(K[0,0],np.complex128)
#=======================================================================

def test_negative_time_fresnel():
#=======================================================================
    "Check that kernel is constructed properly for negative times using the fresnel method"

    tneg = np.linspace(-5,0,50)
    tpos = np.linspace(0,5,50) 
    r = np.linspace(2,6,70)
    Kneg = dipolarkernel(tneg,r,method='fresnel')
    Kpos = dipolarkernel(tpos,r,method='fresnel')

    assert np.all(abs(Kneg - np.flipud(Kpos))<1e-12)
#=======================================================================


def test_negative_time_grid():
#=======================================================================
    "Check that kernel is constructed properly for negative times using the grid method"

    tneg = np.linspace(-5,0,50)
    tpos = np.linspace(0,5,50) 
    r = np.linspace(2,6,70)
    Kneg = dipolarkernel(tneg,r,method='grid')
    Kpos = dipolarkernel(tpos,r,method='grid')

    delta = abs(Kneg - np.flipud(Kpos))

    assert np.all(delta<1e-12)
#=======================================================================


def test_negative_time_integral():
#=======================================================================
    "Check that kernel is constructed properly for negative times using the integral method"

    tneg = np.linspace(-5,0,50)
    tpos = np.linspace(0,5,50) 
    r = 5
    Kneg = dipolarkernel(tneg,r,method='integral')
    Kpos = dipolarkernel(tpos,r,method='integral')

    delta = abs(Kneg - np.flipud(Kpos))

    assert np.allclose(Kneg,np.flipud(Kpos))
#=======================================================================

def test_value_fresnel():
#=======================================================================
    "Test whether kernel matrix element (calculated using Fresnel integrals) is correct."

    # Generate kernel numerically
    t = 1 # µs
    r = 1 # nm
    K = dipolarkernel(t,r,method='fresnel')

    # Kernel value for 1µs and 1nm computed using Mathematica (FresnelC and FresnelS) 
    # and CODATA 2018 values for ge, µB, µ0, and h.
    Kref = 0.024697819895260188

    assert abs(K-Kref) < 1e-14
#=======================================================================

def test_complex_value_fresnel():
#=======================================================================
    "Test whether comple-valued kernel matrix element (calculated using Fresnel integrals) is correct."

    # Generate kernel numerically
    t = 1 # µs
    r = 1 # nm
    K = dipolarkernel(t,r,method='fresnel',complex=True)

    # Kernel value for 1µs and 1nm computed using Mathematica (FresnelC and FresnelS) 
    # and CODATA 2018 values for ge, µB, µ0, and h.
    Kref = 0.0246978198952619 + 1j*0.01380433055548875

    assert abs(K-Kref) < 1e-14
#=======================================================================


def test_complex_value_fresnel_negative_time():
#=======================================================================
    "Test whether comple-valued kernel matrix element (calculated using Fresnel integrals) is correct."

    # Generate kernel numerically
    t = -1 # µs
    r = 1 # nm
    K = dipolarkernel(t,r,method='fresnel',complex=True)

    # Kernel value for 1us and 1nm computed using Mathematica (FresnelC and FresnelS) 
    # and CODATA 2018 values for ge, muB, mu0, and h.
    Kref = 0.0246978198952619 - 1j*0.01380433055548875

    assert abs(K-Kref) < 1e-14
#=======================================================================

def test_value_grid():
#=======================================================================
    "Test whether kernel matrix element (calculated using the grid method) is correct."

    # Generate kernel numerically
    t = 1 # µs
    r = 1 # nm
    K = dipolarkernel(t,r,method='grid',gridsize=2e6)

    # Kernel value for 1us and 1nm computed using Mathematica (FresnelC and FresnelS) 
    # and CODATA 2018 values for ge, muB, mu0, and h.
    Kref = 0.024697819895260188
    print(K.shape)
    assert abs(K-Kref) < 1e-6
#=======================================================================


def test_complex_value_grid():
#=======================================================================
    "Test whether comple-valued kernel matrix element (calculated using the grid method) is correct."

    # Generate kernel numerically
    t = 1 # µs
    r = 1 # nm
    K = dipolarkernel(t,r,method='grid',gridsize=2e6,complex=True)

    # Kernel value for 1us and 1nm computed using Mathematica (FresnelC and FresnelS) 
    # and CODATA 2018 values for ge, muB, mu0, and h.
    Kref = 0.0246978198952619 + 1j*0.01380433055548875

    assert abs(K-Kref) < 1e-6
#=======================================================================


def test_complex_value_grid_negative_time():
#=======================================================================
    "Test whether comple-valued kernel matrix element (calculated using the grid method) is correct."

    # Generate kernel numerically
    t = -1 # µs
    r = 1 # nm
    K = dipolarkernel(t,r,method='grid',gridsize=2e6,complex=True)

    # Kernel value for 1us and 1nm computed using Mathematica (FresnelC and FresnelS) 
    # and CODATA 2018 values for ge, muB, mu0, and h.
    Kref = 0.0246978198952619 - 1j*0.01380433055548875

    assert abs(K-Kref) < 1e-6
#=======================================================================

def test_value_integral():
#=======================================================================
    "Test whether kernel matrix element (calculated using the integal method) is correct."

    # Generate kernel numerically
    t = 1 # µs
    r = 1 # nm
    K = dipolarkernel(t,r,method='integral')

    # Kernel value for 1us and 1nm computed using Mathematica (FresnelC and FresnelS) 
    # and CODATA 2018 values for ge, muB, mu0, and h.
    Kref = 0.024697819895260188

    assert abs(K-Kref) < 1e-7
#=======================================================================

def test_lambda():
#=======================================================================
    "Check that dipolar kernel with modulation depth works"

    t = np.linspace(0,4,150) # µs
    r = np.linspace(1,5,150) # nm
    lam = 0.4
    Kref = ((1-lam) + lam*dipolarkernel(t,r, integralop=False))
    K = dipolarkernel(t,r,mod=lam, integralop=False)

    assert np.allclose(K,Kref)
#=======================================================================

def test_background():
#=======================================================================
    "Check that dipolar kernel builds with the background included"
    t = np.linspace(0,3,80)
    r = np.linspace(2,6,100)
    B = bg_exp(t,0.5)

    lam = 0.25
    Kref = (1 - lam + lam*dipolarkernel(t,r, integralop=False))*B[:,np.newaxis]
    K = dipolarkernel(t,r,mod=lam, bg=B, integralop=False)
    
    assert np.allclose(K,Kref)
#=======================================================================

def test_nobackground():
#=======================================================================
    "Check that dipolar kernel builds correctly without a background specified"
    t = np.linspace(0,3,80)
    r = np.linspace(2,6,100)

    Kref = dipolarkernel(t,r)
    K = dipolarkernel(t,r,bg=None)

    assert np.allclose(K,Kref)
#=======================================================================

def test_multipath():
#=======================================================================
    "Check that multi-pathway kernels are properly generated"

    r = np.linspace(2,6,50)
    t1 = np.linspace(0,10,300)
    t2 = 0.3
    tau1 = 4.24
    tau2 = 4.92
    t = (tau1 + tau2) - (t1 + t2)
    prob = 0.8
    lam = [prob**2, prob*(1-prob)]
    T0 = [0, tau2-t2]

    paths = []
    paths.append({'amp':1-prob})
    paths.append({'amp':prob**2, 'reftime': 0})
    paths.append({'amp':prob*(1-prob), 'reftime': tau2-t2})

    K = dipolarkernel(t,r,pathways=paths, integralop=False)

    Kref = 1-prob
    for p in range(len(lam)):
        Kref = Kref + lam[p]*elementarykernel_twospin(t-T0[p],r,'fresnel',[],[],[ge,ge],None,False)

    assert np.allclose(K,Kref)
#=======================================================================

def test_multipath_background():
#=======================================================================
    "Check that multi-pathway kernels are properly generated with background"

    r = np.linspace(2,6,50)
    t1 = np.linspace(0,10,300)
    t2 = 0.3
    tau1 = 4.24
    tau2 = 4.92
    t = (tau1 + tau2) - (t1 + t2)
    prob = 0.8
    lam = [prob**2, prob*(1-prob)]
    T0 = [0, tau2-t2]
    conc = 50
    Bmodel = lambda t,lam: bg_hom3d(t,conc,lam)

    # Reference
    Kref = 1-prob
    for p in range(len(lam)):
        Kref = Kref + lam[p]*elementarykernel_twospin(t-T0[p],r,'fresnel',[],[],[ge,ge],None,False)
    Bref = 1
    for p in range(len(lam)):
        Bref = Bref*Bmodel((t-T0[p]),lam[p])
    Kref = Kref*Bref[:,np.newaxis]

    paths = []
    paths.append({'amp': 1-prob})
    paths.append({'amp': prob**2, 'reftime': 0})
    paths.append({'amp': prob*(1-prob), 'reftime': tau2-t2})

    # Output
    K = dipolarkernel(t,r,pathways=paths, bg=Bmodel, integralop=False)

    assert np.allclose(K,Kref)
#=======================================================================


def test_multipath_harmonics():
#=======================================================================
    "Check that multi-pathway kernels with higher harmonics are properly generated"

    r = np.linspace(2,6,50)
    t1 = np.linspace(0,10,300)
    t2 = 0.3
    tau1 = 4.24
    tau2 = 4.92
    t = (tau1 + tau2) - (t1 + t2)
    prob = 0.8
    lam = [prob**2, prob*(1-prob)]
    tref = [0, tau2-t2]
    n = [2, 3]

    paths = []
    paths.append({'amp': 1-prob})
    paths.append({'amp': lam[0], 'reftime': tref[0], 'harmonic': n[0]})
    paths.append({'amp': lam[1], 'reftime': tref[1], 'harmonic': n[1]})

    K = dipolarkernel(t,r,pathways=paths, integralop=False)

    Kref = 1-prob
    for p in range(len(lam)):
        Kref = Kref + lam[p]*elementarykernel_twospin(n[p]*(t - tref[p]),r,'fresnel',[],[],[ge,ge],None,False)

    assert np.allclose(K,Kref)
#=======================================================================

def test_excbandwidth_inf_grid():
#=======================================================================
    "Check that specifying a bandwidth effectively changes the kernel (grid method)"

    r = 2.5 # nm
    excitewidth = 15 # MHz
    t = np.linspace(0,2,501) # µs

    # Reference kernel with infinite bandwidth
    Kinf = dipolarkernel(t,r,excbandwidth=inf,method='grid')

    # Kernel with limited bandwidth
    Klim = dipolarkernel(t,r,excbandwidth=excitewidth,method='grid')

    assert not np.array_equal(Kinf,Klim)
#=======================================================================

def test_excbandwidth_inf_integral():
#=======================================================================
    "Check that specifying a bandwidth effectively changes the kernel (integral method)"

    r = 2.5 # nm
    excitewidth = 15 # MHz
    t = np.linspace(0,2,501) # µs

    # Reference kernel with infinite bandwidth
    Kinf = dipolarkernel(t,r,excbandwidth=inf,method='integral')

    # Kernel with limited bandwidth
    Klim = dipolarkernel(t,r,excbandwidth=excitewidth,method='integral')

    assert not np.array_equal(Kinf,Klim)
#=======================================================================

def test_excbandwidth():
#=======================================================================
    """Check kernel with limited excitation bandwidth is numerically correct"""

    r = 2.5 # nm
    excitewidth = 15 # MHz
    t = np.linspace(0,2,501) # µs

    # Kernel using numerical integrals and approx. bandwidth treatment
    K = dipolarkernel(t,r,excbandwidth=excitewidth,method='grid',gridsize=5000)

    # Manual calculation using numerical powder average (reference)
    gridsize = 5001
    costheta = np.linspace(0,1,gridsize)
    wdd = 2*pi*52.04/r**3 # Mrad/s
    q = 1-3*costheta**2
    D_ = 0
    for orientation in q:
        w = wdd*orientation
        D_ = D_ + np.cos(w*abs(t))*np.exp(-w**2/excitewidth**2)
    K0 = D_/gridsize
    K0 = K0.reshape((len(t),1))

    assert np.max(K0-K)<1e-4
#=======================================================================

def test_gvalues():
#=======================================================================
    """Check whether K matrix elements scale properly g-values"""

    t = 2 # nm
    r = 3 # µs
    # isotropic g values of two spins
    g1 = 2.1
    g2 = 2.4

    ge = 2.00231930436256 # free-electron g factor (CODATA 2018 value)

    K1 = dipolarkernel(t,r,g=[g1, g2])
    K2 = dipolarkernel(t,r/(g1*g2/ge**2)**(1/3))

    assert np.max(K1 - K2) < 1e-15
#=======================================================================

def test_r_scaling():
#=======================================================================
    """Check whether K matrix elements scale properly with t and r"""

    t = 2 # nm
    r = 3 # µs
    c = 1.2 # distance scaling factor

    K1 = dipolarkernel(t,r)
    K2 = dipolarkernel(t*c**3,r*c)

    assert np.max(K1 - K2) < 1e-15
#=======================================================================


def test_integralop():
#=======================================================================
    """Check that integral operator-nature of kernel can be disabled"""

    t = np.linspace(0,5,101)
    r = np.linspace(2,5,101)
    dr = np.ones(len(r))*np.mean(np.diff(r))
    dr[0] = dr[0]/2
    dr[-1] = dr[-1]/2
    Kref = dipolarkernel(t,r, integralop=True)/dr
    K = dipolarkernel(t,r, integralop=False)

    assert np.max(K - Kref) < 1e-14
#=======================================================================


def test_nonuniform_r():
#=======================================================================
    "Check that normalization is correct when using a non-uniform distance axis"

    t = 0
    r = np.sqrt(np.linspace(1,7**2,200))
    P = dd_gauss(r,3,0.5)
    K = dipolarkernel(t,r)
    V0 = K@P
    assert np.round(V0,3) == 1
#=======================================================================

def test_orisel_uni_grid():
#=======================================================================
    "Check that orientation selection works for a uniform distribution"

    t = 1
    r = 1
    Ptheta = lambda theta: np.ones_like(theta)

    Kref = dipolarkernel(t,r,method='grid',gridsize=1e5)
    K = dipolarkernel(t,r,method='grid',gridsize=1e5,orisel=Ptheta)
    
    assert np.max(K - Kref) < 1e-10
#=======================================================================

def test_orisel_uni_integral():
#=======================================================================
    "Check that orientation selection works for a uniform distribution"

    t = 1
    r = 1
    Ptheta = lambda theta: np.ones_like(theta)

    Kref = dipolarkernel(t,r,method='integral')
    K = dipolarkernel(t,r,method='integral',orisel=Ptheta)
    
    assert np.max(K - Kref) < 1e-10
#=======================================================================

def test_orisel_value_grid():
#=======================================================================
    "Check that orientation selection works for a uniform distribution"

    t = 1
    r = 1
    thetamean = pi/4
    sigma = pi/3
    Ptheta = lambda theta: 1/sigma/np.sqrt(2*pi)*np.exp(-(theta-thetamean)**2/2/sigma**2)

    # Kernel value for 1us and 1nm computed using Mathematica
    # and CODATA 2018 values for ge, muB, mu0, and h
    Kref = 0.02031864642707554

    K = dipolarkernel(t,r,method='grid',gridsize=1e7,orisel=Ptheta)

    assert abs(K - Kref) < 1e-4
#=======================================================================

def test_orisel_value_integral():
#=======================================================================
    "Check that orientation selection works for a uniform distribution"

    t = 1
    r = 1
    thetamean = pi/4
    sigma = pi/3
    Ptheta = lambda theta: 1/sigma/np.sqrt(2*pi)*np.exp(-(theta-thetamean)**2/2/sigma**2)

    # Kernel value for 1us and 1nm computed using Mathematica
    # and CODATA 2018 values for ge, muB, mu0, and h
    Kref = 0.02031864642707554

    K = dipolarkernel(t,r,method='integral',orisel=Ptheta)
    
    assert abs(K - Kref) < 1e-4
#=======================================================================

def test_docstring():
# ======================================================================
    "Check that the docstring includes all variables and keywords."
    assert_docstring(dipolarkernel)
# ======================================================================

import pytest 

def test_memory_limit():
# ======================================================================
    "Check that the memory limiter works"
    # Construct 12GB kernel (should fit in RAM)
    t = np.linspace(0,5,int(3e4))
    r = np.linspace(1,6,int(5e4))
    with pytest.raises(MemoryError):
        K = dipolarkernel(t,r)
# ======================================================================

def test_twodimensional_size():
#=======================================================================
    "Check the contrusction of kernels with two time dimensions"

    t1 = np.linspace(0,5,10)
    t2 = np.linspace(0,3,20)
    r = np.array(3.5)
    pathways= [
    {'amp': 0.7},
    {'reftime': [0,0], 'amp': 0.3, 'harmonic': [1,0]}    
    ]
    K = dipolarkernel([t1,t2],r,pathways=pathways, method='grid', gridsize=5)

    assert np.allclose(np.shape(K), np.array([10,20,1]))
#=======================================================================

def test_twodimensional_twospins_t1modulated():
#=======================================================================
    "Check the contrusction of kernels with two time dimensions on two spin systems"

    t1 = np.linspace(0,5,50)
    t2 = np.linspace(0,3,50)
    r = np.array(3.5)
    pathways= [
        {'amp': 0.7},
        {'reftime': [0,0], 'amp': 0.3, 'harmonic': [1,0]}    
    ]
    K = dipolarkernel([t1,t2],r,pathways=pathways)
    Kref = dipolarkernel(t1,r,mod=0.3)

    assert np.allclose(K.squeeze(),Kref)
#=======================================================================

def test_twodimensional_twospins_t2modulated():
#=======================================================================
    "Check the contrusction of kernels with two time dimensions on two spin systems"

    t1 = np.linspace(0,5,50)
    t2 = np.linspace(0,3,50)
    r = np.array(3.5)
    pathways= [
        {'amp': 0.7},
        {'reftime': [0,0], 'amp': 0.3, 'harmonic': [0,1]}    
    ]
    K = dipolarkernel([t1,t2],r,pathways=pathways)
    Kref = dipolarkernel(t2,r,mod=0.3)

    assert np.allclose(K.squeeze().T,Kref)
#=======================================================================

def test_threespin_matrixshape():
#=======================================================================
    "Check that three-spin kernels are still bidimensional"

    t = np.linspace(0,5,20)
    r1 = np.linspace(1,4,10)
    r2 = np.linspace(1,4,10)
    r3 = np.linspace(1,4,10)
    pathways= [
    {'amp': 0.7},
    {'reftime': (0,0,None),'amp': 0.3, 'harmonic': (1,1,0)}    
    ]
    K = dipolarkernel(t,[r1,r2,r3],pathways=pathways)

    assert np.allclose(np.shape(K), np.array([20,10]))
#=======================================================================

def test_threespin_grid_value():
#=======================================================================
    "Check that the three-spin kernels are accurate"

    t = np.array(0.1)
    r1 = np.array(3.5)
    r2 = np.array(3.2)
    r3 = np.array(4.2)
    pathways= [{'reftime': (0,0,None),'amp': 1, 'harmonic': (1,1,0)}]
    K = dipolarkernel(t,[r1,r2,r3],pathways=pathways,gridsize=15000)

    Kref = 0.5136853637 # Computed with Mathematica 11.3 (NIntegrate)

    assert abs(K - Kref)<1e-2
#=======================================================================

def test_threespin_grid_complexvalue():
#=======================================================================
    "Check that the three-spin kernels are accurate"

    t = np.array(0.25)
    r1 = np.array(3.5)
    r2 = np.array(3.2)
    r3 = np.array(4.2)
    pathways= [{'reftime': (0,0,None),'amp': 1, 'harmonic': (1,1,0)}]
    K = dipolarkernel(t,[r1,r2,r3],pathways=pathways,gridsize=15000, complex=True)

    Kref = -0.0236232 + 0.109215j  # Computed with Mathematica 11.3 (NIntegrate)

    assert abs(K - Kref)<1e-2
#=======================================================================

def test_threespin_twodimensional():
#=======================================================================
    "Check that three-spin and bidimensional kernels are of the correct size"

    t1 = np.linspace(0,5,20)
    t2 = np.linspace(0,5,30)
    r1 = np.linspace(1,4,10)
    r2 = np.linspace(1,4,10)
    r3 = np.linspace(1,4,10)
    pathways= [
    {'amp': 0.7},
    {'reftime': ([0,0],[0,0],[None,None]),'amp': 0.3, 'harmonic': ([1,0],[0,1],[0,0])}    
    ]
    K = dipolarkernel([t1,t2],[r1,r2,r3],pathways=pathways)

    assert np.allclose(np.shape(K), np.array([20,30,10]))
#=======================================================================

def test_args_arraytype_mod():
#=======================================================================
    "Check that input arguments are insensitive to type of array"
    t = 0.5
    r = 5
    mod = 0.3
    K1 = dipolarkernel(t,r,mod=mod)
    mod = np.array(0.3)
    K2 = dipolarkernel(t,r,mod=mod)
    mod = np.array([0.3])
    K3 = dipolarkernel(t,r,mod=mod)

    assert K1==K2 and K2==K3
#=======================================================================

def test_args_arraytype_t():
#=======================================================================
    "Check that input arguments are insensitive to type of array"
    r = 5
    t = 0.3
    K1 = dipolarkernel(t,r)
    t = np.array(0.3)
    K2 = dipolarkernel(t,r)
    t = np.array([0.3])
    K3 = dipolarkernel(t,r)

    assert K1==K2 and K2==K3
#=======================================================================

def test_args_arraytype_r():
#=======================================================================
    "Check that input arguments are insensitive to type of array"
    t = 0.3
    r = 5
    K1 = dipolarkernel(t,r)
    r = np.array(5)
    K2 = dipolarkernel(t,r)
    r = np.array([5])
    K3 = dipolarkernel(t,r)

    assert K1==K2 and K2==K3
#=======================================================================

def test_interpolation_value():
#=======================================================================
    "Test whether kernel matrix element (calculated by interpolation speedup) is correct."

    # Generate kernel numerically
    t = 1 # µs
    r = 1 # nm
    tinterp = np.arange(0,2,0.005) # µs
    K = dipolarkernel(t,r,method='fresnel',tinterp=tinterp)

    # Kernel value for 1µs and 1nm computed using Mathematica (FresnelC and FresnelS) 
    # and CODATA 2018 values for ge, µB, µ0, and h.
    Kref = 0.024697819895260188

    assert abs(K-Kref) < 1e-14
#=======================================================================

def test_interpolation_matrix_shape():
#=======================================================================
    "Test whether interpolation returns the correct kernel shape"

    t = np.linspace(0,5,40) # µs
    r = np.linspace(1,4,30) # nm
    tinterp = t # µs
    Kref = dipolarkernel(t,r,method='fresnel')
    K = dipolarkernel(t,r,method='fresnel',tinterp=tinterp)

    assert np.allclose(K.shape,Kref.shape)
#=======================================================================

def test_interpolation_vector():
#=======================================================================
    "Test whether interpolation returns the same vector as without"

    t = np.linspace(0,5,100) # µs
    r = 4 # nm
    tinterp = t # µs
    Kref = dipolarkernel(t,r,method='fresnel')
    K = dipolarkernel(t,r,method='fresnel',tinterp=tinterp)

    assert np.allclose(K,Kref)
#=======================================================================

def test_interpolation_reftimes():
#=======================================================================
    "Test whether interpolation returns the corrected time-shifted vector"

    t = np.linspace(0,5,100) # µs
    r = 4 # nm
    tinterp = np.linspace(-2,7,400) # µs
    pathways = [{'amp':1, 'reftime':1}]
    Kref = dipolarkernel(t,r,pathways=pathways, method='fresnel')
    K = dipolarkernel(t,r,pathways=pathways, method='fresnel',tinterp=tinterp)

    assert np.allclose(K,Kref,rtol=1e-5)
#=======================================================================

def test_interpolation_threespin():
#=======================================================================
    "Check that the three-spin kernels are accurate"

    t = np.linspace(0,5,100)
    r1 = np.array(3.5)
    r2 = np.array(3.2)
    r3 = np.array(4.2)
    tinterp = np.linspace(-1,4,500)
    pathways= [{'reftime': (None,1,None),'amp': 1, 'harmonic': (0,1,0)}]
    Kref = dipolarkernel(t,[r1,r2,r3],pathways=pathways)
    K = dipolarkernel(t,[r1,r2,r3],pathways=pathways,tinterp=tinterp)

    assert np.allclose(K,Kref,rtol=1e-5)
#=======================================================================

def test_interpolation_outofrange():
#=======================================================================
    "Test whether interpolation rolls back if outside of range"

    t = np.linspace(0,5,40) # µs
    r = np.linspace(1,4,30) # nm
    tinterp = np.linspace(0,4,40) # µs
    Kref = dipolarkernel(t,r,method='fresnel')
    K = dipolarkernel(t,r,method='fresnel',tinterp=tinterp)

    assert np.allclose(K.shape,Kref.shape)
#=======================================================================
