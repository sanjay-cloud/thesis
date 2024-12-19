import numpy as np

def G03(x):
    n = len(x)

    obj = (np.sqrt(n) ** n) * np.prod(x)
    
    g1 = np.sum(np.square(x)) - 1
    
    return np.array([obj, g1])