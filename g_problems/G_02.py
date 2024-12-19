import numpy as np

def G02(x):
    # Objective function
    obj = np.abs((np.sum(np.cos(x)**4) - 2 * np.prod(np.cos(x)**2))/  np.sqrt(np.sum(np.arange(1,len(x)+1) * x**2)))
    
    # Constraints
    g1 = 0.75 - np.prod(x)
    g2 = np.sum(x) - 7.5 * len(x)
    
    return np.array([-obj, g1, g2])
