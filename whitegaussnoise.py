import numpy as np
import matrixtools as mat

def whitegaussnoise(t,level):
    N = len(t)

    noise = np.random.normal(0,level,N)
    noise = mat.column(noise)
    return noise
