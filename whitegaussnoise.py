import numpy as np

def whitegaussnoise(t,level):
    N = len(t)

    noise = np.random.normal(0,level,N)
    return noise
