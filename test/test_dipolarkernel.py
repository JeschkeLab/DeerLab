
import numpy as np
from numpy import pi, inf, NaN
from deerlab.bg_models import bg_hom3d,bg_exp
from deerlab.dd_models import dd_gauss
from deerlab.dipolarkernel import dipolarkernel,elementarykernel
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

    delta = abs(Kneg - np.flipud(Kpos))

    assert np.all(delta<1e-12)
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
    r = np.linspace(2,6,10)
    Kneg = dipolarkernel(tneg,r,method='integral')
    Kpos = dipolarkernel(tpos,r,method='integral')

    delta = abs(Kneg - np.flipud(Kpos))

    assert np.all(delta<1e-12)
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
    K = dipolarkernel(t,r,method='grid',nKnots=2e6)

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
    K = dipolarkernel(t,r,method='grid',nKnots=2e6,complex=True)

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
    K = dipolarkernel(t,r,method='grid',nKnots=2e6,complex=True)

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

    assert np.all(abs(K-Kref) < 1e-14)
#=======================================================================

def test_background():
#=======================================================================
    "Check that dipolar kernel builds with the background included"
    t = np.linspace(0,3,80)
    r = np.linspace(2,6,100)
    B = bg_exp(t,0.5)

    lam = 0.25
    KBref = (1 - lam + lam*dipolarkernel(t,r, integralop=False))*B[:,np.newaxis]
    KB = dipolarkernel(t,r,mod=lam, bg=B, integralop=False)
    
    assert np.all(abs(KB - KBref) < 1e-5)
#=======================================================================

def test_nobackground():
#=======================================================================
    "Check that dipolar kernel builds correctly without a background specified"
    t = np.linspace(0,3,80)
    r = np.linspace(2,6,100)

    K1 = dipolarkernel(t,r)
    
    K2 = dipolarkernel(t,r,bg=None)
    assert np.all(abs(K2 - K1) < 1e-5)
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
    paths.append([1-prob])
    paths.append([prob**2, 0])
    paths.append([prob*(1-prob), tau2-t2])

    K = dipolarkernel(t,r,pathways=paths, integralop=False)

    Kref = 1-prob
    for p in range(len(lam)):
        Kref = Kref + lam[p]*elementarykernel(t-T0[p],r,'fresnel',[],[],[ge,ge],None,False)

    assert np.all(abs(K-Kref) < 1e-3)
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
        Kref = Kref + lam[p]*elementarykernel(t-T0[p],r,'fresnel',[],[],[ge,ge],None,False)
    Kref = Kref
    
    Bref = 1
    for p in range(len(lam)):
        Bref = Bref*Bmodel((t-T0[p]),lam[p])
    KBref = Kref*Bref[:,np.newaxis]

    paths = []
    paths.append([1-prob])
    paths.append([prob**2, 0])
    paths.append([prob*(1-prob), tau2-t2])

    # Output
    KB = dipolarkernel(t,r,pathways=paths, bg=Bmodel, integralop=False)

    assert np.all(abs(KB - KBref) < 1e-3)
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
    T0 = [0, tau2-t2]
    n = [2, 3]

    paths = []
    paths.append([1-prob])
    paths.append([prob**2, 0, 2])
    paths.append([prob*(1-prob), tau2-t2,3])

    K = dipolarkernel(t,r,pathways=paths, integralop=False)

    Kref = 1-prob
    for p in range(len(lam)):
        Kref = Kref + lam[p]*elementarykernel(n[p]*(t-T0[p]),r,'fresnel',[],[],[ge,ge],None,False)
    Kref = Kref

    assert np.max(K-Kref) < 1e-3
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
    K = dipolarkernel(t,r,excbandwidth=excitewidth,method='grid')

    # Manual calculation using numerical powder average (reference)
    nKnots = 5001
    costheta = np.linspace(0,1,nKnots)
    wdd = 2*pi*52.04/r**3 # Mrad/s
    q = 1-3*costheta**2
    D_ = 0
    for orientation in q:
        w = wdd*orientation
        D_ = D_ + np.cos(w*abs(t))*np.exp(-w**2/excitewidth**2)
    K0 = D_/nKnots
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

    Kref = dipolarkernel(t,r,method='grid',nKnots=1e5)
    K = dipolarkernel(t,r,method='grid',nKnots=1e5,orisel=Ptheta)
    
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

    K = dipolarkernel(t,r,method='grid',nKnots=1e7,orisel=Ptheta)
    print(K/Kref)
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