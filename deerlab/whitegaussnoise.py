import numpy as np

def whitegaussnoise(t,level=1,rescale=False,seed=None):
    r"""
    Generates a vector of white Gaussian (normal) noise
    
    Parameters
    ----------
    t : array_like
        Time axis.
    level : float scalar
        Noise level, i.e. standard deviation of underlying Gaussian distribution.
    rescale : boolean, optional
        If true, rescales the noise vector such that its standard deviation is exactly equal
        to ``level``. If false (default), the standard deviation of the noise vector can deviate
        slightly from ``level``, particularly for short vectors.
    seed : integer scalar, optional
        If ``None`` (default), do not seed the random number generator. If an integer scalar is
        given (e.g. ``seed=137``), seed the random number generator with this number.
    
    Returns
    -------
    noise : ndarray
        Noise vector.
    
    Notes
    -----
    The noise vector is generated using pseudo-random numbers generated with NumPy. 
    Without seeding the random number generator, subsequent calls of ``whitegaussnoise`` return
    different realizations of the noise vector. To obtain a reproducible noise realization, seed
    the random number generator by using the ``seed`` kewyword arguement, or call ``numpy.random.seed(k)``
    with some integer number ``k`` before calling ``whitegaussnoise``.
    """
    
    # Seed RNG if wanted
    if seed is not None:
        np.random.seed(seed)
    
    # Draw from standard normal distribution
    N = len(np.atleast_1d(t))
    noise = np.random.normal(0,level,N)
    
    # Rescale to sample std if wanted
    if rescale:
        noise = level/np.std(noise)*noise
    
    return noise
