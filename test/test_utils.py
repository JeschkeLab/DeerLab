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
