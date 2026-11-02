import numpy as np
from deerlab import sophegrid


# ======================================================================
# Test sophegrid function

def test_sophegrid():
    # Test that the function returns the values and weights summing to 1
    # Comparison taken from EasySpin 6.0
    phi, theta, weights = sophegrid(4,np.pi*2,3)

    assert np.allclose(weights.sum(),1)
    assert np.allclose(phi, np.array([0,0,1.57079632679490,3.14159265358979,4.71238898038469,0,0.785398163397448,1.57079632679490,2.35619449019235,3.14159265358979,3.92699081698724,4.71238898038469,5.49778714378214]))
    assert np.allclose(theta, np.array([0,0.785398163397448,0.785398163397448,0.785398163397448,0.785398163397448,1.57079632679490,1.57079632679490,1.57079632679490,1.57079632679490,1.57079632679490,1.57079632679490,1.57079632679490,1.57079632679490]))
    assert np.allclose(weights*4*np.pi, np.array([0.956558005801449,1.70021769237074,1.70021769237074,1.70021769237074,1.70021769237074,0.601117729884346,0.601117729884346,0.601117729884346,0.601117729884346,0.601117729884346,0.601117729884346,0.601117729884346,0.601117729884346]))

    phi, theta, weights = sophegrid(1,np.pi/2,5,closed_phi=True)
    
    assert np.allclose(weights.sum(),1)
    assert np.allclose(phi, np.array([0, 0, 1.5708, 0, 0.7854, 1.5708, 0, 0.5236, 1.0472, 1.5708, 0, 0.3927, 0.7854, 1.1781, 1.5708]),rtol=1e-4)
    assert np.allclose(theta, np.array([0, 0.3927, 0.3927, 0.7854, 0.7854, 0.7854, 1.1781, 1.1781, 1.1781, 1.1781, 1.5708, 1.5708, 1.5708, 1.5708, 1.5708]),rtol=1e-4)
    assert np.allclose(weights*4*np.pi, np.array([0.24146, 0.93818, 0.93818, 0.86676, 1.7335, 0.86676, 0.75499, 1.51, 1.51, 0.75499, 0.30645, 0.61289, 0.61289, 0.61289, 0.30645]),rtol=1e-4)


    phi, theta, weights = sophegrid(0,0,6)    
    assert np.allclose(weights.sum(),1)
    assert np.allclose(phi, np.array([0,	0,	0,	0,	0,	0]),rtol=1e-3)
    assert np.allclose(theta, np.array([0,	0.3142,	0.6283,	0.9425,	1.2566,	1.5708]),rtol=1e-3)
    assert np.allclose(weights*4*np.pi, np.array([0.1547,	1.2149,	2.3110,	3.1808,	3.7392,	1.9658]),rtol=1e-3)

    phi, theta, weights = sophegrid(-1,0,6)    

    assert np.allclose(weights.sum(),1)
    assert np.allclose(phi, np.array([0]),rtol=1e-4)
    assert np.allclose(theta, np.array([0]),rtol=1e-4)
    assert np.allclose(weights*4*np.pi, np.array([12.5664]),rtol=1e-4)

def test_choleskycovmat_roundtrip():
    # Test that the factors of a known covariance matrix, ordered as in dipolarmodel
    # (ℓ11,...,ℓQQ, ℓ21,ℓ31,...,ℓQ1, ℓ32,...), give back the same covariance matrix
    from deerlab.utils import choleskycovmat
    for Q in [3,4,6]:
        rng = np.random.default_rng(Q)
        A = rng.normal(size=(Q,Q))
        Σ = 0.1*(A@A.T + Q*np.eye(Q))
        L = np.linalg.cholesky(Σ)
        cholfactors = [L[q,q] for q in range(Q)] + [L[j,q] for q in range(Q) for j in range(q+1,Q)]

        assert np.allclose(choleskycovmat(Q,cholfactors),Σ)

def test_choleskycovmat_named_factors():
    # Test that each named factor ℓjq is placed at the lower triangular position (j,q)
    from deerlab.utils import choleskycovmat
    # Only ℓ21 non-zero: correlation between dimensions 1 and 2 only
    Σ = choleskycovmat(3,[0.3,0.3,0.3, 0.2,0,0])
    assert Σ[0,1]!=0 and Σ[0,2]==0 and Σ[1,2]==0
    # Only ℓ32 non-zero: correlation between dimensions 2 and 3 only
    Σ = choleskycovmat(3,[0.3,0.3,0.3, 0,0,0.2])
    assert Σ[1,2]!=0 and Σ[0,1]==0 and Σ[0,2]==0
    # Docstring example
    assert np.allclose(choleskycovmat(3,[3.0,2.0,1.0,0.5,0.5,0.5]), [[9,1.5,1.5],[1.5,4.25,1.25],[1.5,1.25,1.5]])

def test_config():
    from deerlab.utils import show_config

    config_dict = show_config(mode='dicts')

    assert 'python_version' in config_dict
    assert 'deerlab_version' in config_dict
    assert 'cpu_cores' in config_dict
    assert 'operating_system' in config_dict
    assert 'operating_system_version' in config_dict

    assert 'version' in config_dict['numpy']
    assert 'BLAS' in config_dict['numpy']
    assert 'LAPACK' in config_dict['numpy']
    assert 'version' in config_dict['scipy']
    assert 'BLAS' in config_dict['scipy']
    assert 'LAPACK' in config_dict['scipy']
    
    show_config()   
