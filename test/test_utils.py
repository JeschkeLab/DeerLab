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

    phi, theta, weights = sophegrid(1,np.pi/2,3,closed_phi=True)

    assert np.allclose(weights.sum(),1)
    assert np.allclose(phi, np.array([0, 0, 1.5708, 0, 0.7854, 1.5708]),rtol=1e-4)
    assert np.allclose(theta, np.array([0, 0.7854, 0.7854, 1.5708, 1.5708, 1.5708]),rtol=1e-4)
    assert np.allclose(weights*4*np.pi, np.array([0.9566, 3.4004, 3.4004, 1.2022, 2.4045, 1.2022]),rtol=1e-4)
