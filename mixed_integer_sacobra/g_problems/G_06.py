import numpy as np

def G06(x):
    obj = (x[0] - 10)**3 + (x[1] - 20)**3
    
    g1 = -(x[0] - 5)**2 - (x[1] - 5)**2 + 100
    g2 = (x[0] - 6)**2 + (x[1] - 5)**2 - 82.81
    
    return np.array([obj, g1, g2])
