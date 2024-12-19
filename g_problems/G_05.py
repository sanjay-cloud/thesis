import numpy as np

def G05(x):
    obj =  3*x[0] + 0.000001*x[0]**3 + 2*x[1] + (0.000002/3)*(x[1]**3)
    g1 = -x[3] + x[2] - 0.55
    g2 = -x[2] + x[3] - 0.55

    # Equality constraints
    h3 = 1000*np.sin(-x[2] - 0.25) + 1000*np.sin(-x[3] - 0.25) + 894.8 - x[0]
    h4 = 1000*np.sin(x[2] - 0.25) + 1000*np.sin(x[2] - x[3] - 0.25) + 894.8 - x[1]
    h5 = 1000*np.sin(x[3] - 0.25) + 1000*np.sin(x[3] - x[2] - 0.25) + 1294.8

    return np.array([obj, g1, g2, h3, h4, h5])
