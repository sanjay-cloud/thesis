import numpy as np
from g_problem import GProblem

def problem():
    lower = np.array([100,1000,1000,10,10,10,10,10])
    upper = np.array([10000,10000,10000,1000,1000,1000,1000,1000])
    return GProblem("G10",problem_type_mixint = True, objective = objective, constraints = constraints, 
                    lower_bounds = lower, upper_bounds = upper, dimensions = len(lower), real_indices = [0,3,4,5,6,7],
                    integer_indices=[1,2], minimize = True, xopt=[579.3167, 1359.943,5110.071,182.0174,295.5985,217.9799,286.4162,395.5979],
                    fopt= 7049.330923, x_ticks=2000, discrete_steps=100)

def objective(x):
    x = np.array(x)
    return x[0] + x[1] + x[2]

def constraints(x):
    x = np.array(x)
    g1 = -1 + 0.0025*(x[3] + x[5])
    g2 = -1 + 0.0025*(x[4] + x[6] - x[3])
    g3 = -1 + 0.01*(x[7] - x[4])
    g4 = -x[0]*x[5] + 833.33252*x[3] + 100*x[0] - 83333.333
    g5 = -x[1]*x[6] + 1250*x[4] + x[1]*x[3] - 1250*x[3]
    g6 = -x[2]*x[7] + 1250000 + x[2]*x[4] - 2500*x[4]

    return [g1, g2, g3, g4, g5, g6]


