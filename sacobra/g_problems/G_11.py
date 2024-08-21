import numpy as np

def G11(x):
    obj = (x[0])**2 + (x[1]-1)**2
    
    g1 = x[0]**2 - x[1]
    
    return np.array([obj, g1])