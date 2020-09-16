
import numpy as np
from deerlab import correctzerotime

def test_correction():
#=======================================================================
    "Check that the time-axis is corrected apropiately"

    t_truth = np.arange(-5,80,0.5)
    y = 3000 - t_truth**2
    y = y/max(y)
    t = t_truth + abs(min(t_truth))

    t_corr = correctzerotime(y,t)

    assert max(abs(t_corr - t_truth)) < 1e-10
#=======================================================================

def test_late_times():
#=======================================================================
    "Check that zero-time correction works when maximum is at later times"

    
    t_truth = np.linspace(-5,1,400)
    y = 10 - t_truth**2
    y = y/max(y)
    tshift = 1.2343

    t = t_truth + tshift

    t_corr = correctzerotime(y,t)

    assert max(abs(t_corr - t_truth)) < 1e-10
#=======================================================================

def test_first_element():
#=======================================================================
    "Check that zero-time correction works when maximum is at later times"

    
    t_truth = np.linspace(0,5,400)
    y = 10 - t_truth**2
    y = y/max(y)
    tshift = 1.2343

    t = t_truth + tshift

    t_corr = correctzerotime(y,t)

    assert max(abs(t_corr - t_truth)) < 1e-10
#=======================================================================

def test_last_element():
#=======================================================================
    "Check that zero-time correction works when maximum is at later times"

    
    t_truth = np.linspace(-5,0,400)
    y = 10 - t_truth**2
    y = y/max(y)
    tshift = 1.2343

    t = t_truth + tshift

    t_corr = correctzerotime(y,t)

    assert max(abs(t_corr - t_truth)) < 1e-10
#=======================================================================
