import numpy as np

def G03(x):
    obj = (x[0] - x[1])**2 + (x[1] + x[2] - 2)**2 + (x[3] - 1)**2 + (x[4] - 1)**2
    
    g1 = x[0] + x[1]**2 + x[2]**2 - 1
    g2 = x[0]**2 + x[1] + x[2]**2 - 1
    g3 = x[0]**2 + x[1]**2 + x[2] - 1
    
    return np.array([obj, g1, g2, g3])