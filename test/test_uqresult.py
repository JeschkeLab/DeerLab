import numpy as np 
from deerlab import UQResult, dd_gauss
from deerlab.utils import ovl
from scipy.stats import chi2

# Gaussian uncertainty distribution parameters
std = [0.3,0.6]
means = [5.0,4.0]

# Covariance matrix
covmat = np.diag(np.array(std)**2)

# Bootstrapping
samples = []
for mean,sigma in zip(means,std): 
    np.random.seed(0)
    samples.append(sigma*np.random.randn(200000) + mean)

# Reference percentiles
p95,p50,p5 = np.zeros(2),np.zeros(2),np.zeros(2)
for n,sample in enumerate(samples):
    p95[n] = np.percentile(samples[n],95)
    p5[n] = np.percentile(samples[n],5)
    p50[n] = np.percentile(samples[n],50)

# Reference confidence intervals
ci95,ci90,ci50 = np.zeros((2,2)),np.zeros((2,2)),np.zeros((2,2))
for n,sample in enumerate(samples):
    ci95[n,:] = np.array([np.percentile(samples[n],2.5),np.percentile(samples[n],97.5)])
    ci90[n,:] = np.array([np.percentile(samples[n],5.0),np.percentile(samples[n],95.0)])
    ci50[n,:] = np.array([np.percentile(samples[n],25),np.percentile(samples[n],75)])

# Profile likelihood simulation
x = np.linspace(0,10,1000)
pdf1 = dd_gauss(x,means[0],std[0])
pdf2 = dd_gauss(x,means[1],std[1])
pdf1 /= max(pdf1)
pdf2 /= max(pdf2)
σ = 0.01
obj2likelihood = lambda f: 1/np.sqrt(σ*2*np.pi)*np.exp(-1/2*f/σ**2)
likelihood2obj = lambda L: -2*np.log(L*np.sqrt(σ*2*np.pi))*σ**2
threshold = lambda coverage: σ**2*chi2.ppf(coverage, df=1) + likelihood2obj(max(pdf1))
profile1 = {'y': likelihood2obj(pdf1), 'x':x}
profile2 = {'y': likelihood2obj(pdf2), 'x':x}

# Construct uncertainty quantification objects
uq_covariance = UQResult('covariance',data=np.array(means),covmat=covmat)
uq_bootstrap = UQResult('bootstrap',data=np.vstack(samples).T)
uq_profile = UQResult('profile',data=np.array(means),profiles=[profile1,profile2],threshold=threshold,noiselvl=σ)


# ------------------------------------------------------------
def assert_uq(uq,attr):
    if attr=='mean':
        assert np.allclose(uq.mean,means,rtol=1e-2)
    elif attr=='std':
        assert np.allclose(uq.std,std,rtol=1e-2)
    elif attr=='median':
        assert np.allclose(uq.median,p50,rtol=1e-2)

    elif attr=='ci':
        assert np.allclose(uq.ci(95),ci95,rtol=1e-2)
        assert np.allclose(uq.ci(90),ci90,rtol=1e-2)
        assert np.allclose(uq.ci(50),ci50,rtol=1e-2)

    elif attr=='percentile':
        assert np.allclose(uq.percentile(95),p95,rtol=1e-2)
        assert np.allclose(uq.percentile(5),p5,rtol=1e-2)
        assert np.allclose(uq.percentile(50),p50,rtol=1e-2)

    elif attr=='pardist':
        x1,pdf1 = uq.pardist(0)
        x2,pdf2 = uq.pardist(1)
        xs = [x1,x2]
        pdfs = [pdf1,pdf2]
        pdfs_ref = [dd_gauss(x,mean,sigma) for x,mean,sigma in zip(xs,means,std)]
        assert ovl(pdfs[0],pdfs_ref[0])>0.99
        assert ovl(pdfs[1],pdfs_ref[1])>0.99
# ------------------------------------------------------------

#------------------------------------------------------------
#                          COVARIANCE
#------------------------------------------------------------

# ------------------------------------------------------------
def test_uq_covariance_mean():
    'Check the means of the uncertainty distributions'
    assert_uq(uq_covariance,'mean')
# ------------------------------------------------------------

# ------------------------------------------------------------
def test_uq_covariance_median():
    'Check the medians of the uncertainty distributions'
    assert_uq(uq_covariance,'median')
# ------------------------------------------------------------

# ------------------------------------------------------------
def test_uq_covariance_std():
    'Check the standard deviations of the uncertainty distributions'
    assert_uq(uq_covariance,'std')
# ------------------------------------------------------------

# ------------------------------------------------------------
def test_uq_covariance_ci():
    'Check the confidence intervals of the uncertainty distributions'
    assert_uq(uq_covariance,'ci')
# ------------------------------------------------------------

# ------------------------------------------------------------
def test_uq_covariance_percentile():
    'Check the percentiles of the uncertainty distributions'
    assert_uq(uq_covariance,'percentile')
# ------------------------------------------------------------

# ------------------------------------------------------------
def test_uq_covariance_pardist():
    'Check the estimated uncertainty distributions'
    assert_uq(uq_covariance,'pardist')
# ------------------------------------------------------------

#------------------------------------------------------------
#                          BOOTSTRAP
#------------------------------------------------------------

# ------------------------------------------------------------
def test_uq_bootstrap_mean():
    'Check the means of the uncertainty distributions'
    assert_uq(uq_bootstrap,'mean')
# ------------------------------------------------------------

# ------------------------------------------------------------
def test_uq_bootstrap_median():
    'Check the medians of the uncertainty distributions'
    assert_uq(uq_bootstrap,'median')
# ------------------------------------------------------------

# ------------------------------------------------------------
def test_uq_bootstrap_std():
    'Check the standard deviations of the uncertainty distributions'
    assert_uq(uq_bootstrap,'std')
# ------------------------------------------------------------

# ------------------------------------------------------------
def test_uq_bootstrap_ci():
    'Check the confidence intervals of the uncertainty distributions'
    assert_uq(uq_bootstrap,'ci')
# ------------------------------------------------------------

# ------------------------------------------------------------
def test_uq_bootstrap_percentile():
    'Check the percentiles of the uncertainty distributions'
    assert_uq(uq_bootstrap,'percentile')
# ------------------------------------------------------------

# ------------------------------------------------------------
def test_uq_bootstrap_pardist():
    'Check the estimated uncertainty distributions'
    assert_uq(uq_bootstrap,'pardist')
# ------------------------------------------------------------


#------------------------------------------------------------
#                          PROFILE
#------------------------------------------------------------

# ------------------------------------------------------------
def test_uq_profile_mean():
    'Check the means of the uncertainty distributions'
    assert_uq(uq_profile,'mean')
# ------------------------------------------------------------

# ------------------------------------------------------------
def test_uq_profile_median():
    'Check the medians of the uncertainty distributions'
    assert_uq(uq_profile,'median')
# ------------------------------------------------------------

# ------------------------------------------------------------
def test_uq_profile_std():
    'Check the standard deviations of the uncertainty distributions'
    assert_uq(uq_profile,'std')
# ------------------------------------------------------------

# ------------------------------------------------------------
def test_uq_profile_ci():
    'Check the confidence intervals of the uncertainty distributions'
    assert_uq(uq_profile,'ci')
# ------------------------------------------------------------

# ------------------------------------------------------------
def test_uq_profile_percentile():
    'Check the percentiles of the uncertainty distributions'
    assert_uq(uq_profile,'percentile')
# ------------------------------------------------------------

# ------------------------------------------------------------
def test_uq_profile_pardist():
    'Check the estimated uncertainty distributions'
    assert_uq(uq_profile,'pardist')
# ------------------------------------------------------------
