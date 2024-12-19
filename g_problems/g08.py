import numpy as np
from g_problem import GProblem

def problem():
    lower = np.array([0,0])
    upper = np.array([10,10])
    integer_indices=[0,1]
    real_indices =  np.setdiff1d(np.arange(len(lower)), integer_indices)
    return GProblem("G08", problem_type_mixint = False, objective = objective, constraints = constraints, 
                    lower_bounds = lower, upper_bounds = upper, dimensions = len(lower), real_indices= [1],
                    integer_indices=[0,1], minimize = False, xopt = [1.06875651, 4.26202387], fopt = -0.011182)

def objective(x):
    x = np.array(x)
    return - (np.sin(2 * np.pi * x[0])**3 * np.sin(2 * np.pi * x[1])) / (x[0]**3 * (x[0] + x[1]))

def constraints(x):
    x = np.array(x)
    g1 = x[0]**2 - x[1] + 1
    g2 = 1 - x[0] + (x[1] - 4)**2
    return [g1, g2]


