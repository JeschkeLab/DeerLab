import numpy as np

def whitegaussnoise(t,level,rescale=False):
    r"""White Gaussian noise generator

    Parameters
    ----------
    t : array_like
        Time axis.
    level : float scalar
        Noise level (standard deviation of noise vector)
    rescale : boolean
        Enable/disable rescaling of noise vector to enforce the exact noise level.

    Returns
    -------
    noise : ndarray
        Noise vector.

    Notes
    -----
    The noise vector is generated from a vector of pseudo-random numbers generated with numpy. 
    Each call of ``whitegaussnoise`` will return a different realization of the noise vector. To set the output 
    to a fixed noise realization, the random number generator must be fixed by calling ``numpy.random.seed(k)`` 
    before ``whitegaussnoise``, where ``k`` is some integer number.

    """
    t = np.atleast_1d(t)
    N = len(t)
 
    noise = np.random.normal(0,level,N)

    if rescale:
      noise = level*noise/np.std(noise)

    return noise

