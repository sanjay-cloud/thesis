import numpy as np

def G10(x):
    obj = x[0] + x[1] + x[2]
    
    g1 = -1 + 0.0025*(x[3] + x[5])
    g2 = -1 + 0.0025*(x[4] + x[6] - x[3])
    g3 = -1 + 0.01*(x[7] - x[4])
    g4 = -x[0]*x[5] + 833.33252*x[3] + 100*x[0] - 83333.333
    g5 = -x[1]*x[6] + 1250*x[4] + x[1]*x[3] - 1250*x[3]
    g6 = -x[2]*x[7] + 1250000 + x[2]*x[4] - 2500*x[4]

    return np.array([obj, g1, g2, g3, g4, g5, g6])
