import numpy as np

def column(x):
    x = np.reshape(x,(len(x),1))
    return x


def row(x):
    x = np.reshape(x,(1,len(x)))
    return x