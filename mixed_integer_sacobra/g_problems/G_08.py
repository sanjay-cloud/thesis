import numpy as np

def G08(x):
    # Objective function
    obj = - (np.sin(2 * np.pi * x[0])**3 * np.sin(2 * np.pi * x[1])**3) / (x[0]**3 * (x[0] + x[1]))
    
    # Constraint
    g1 = x[0]**2 - x[1] + 1
    g2 = 1 - x[0] + (x[1] - 4)**2
    
    return np.array([obj, g1, g2])