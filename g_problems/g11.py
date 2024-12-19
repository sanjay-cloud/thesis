import numpy as np
from g_problem import GProblem

def problem():
    lower = np.array([-1,-1])
    upper = np.array([1,1])
    return GProblem("G11",problem_type_mixint = True, objective = objective, constraints = constraints, 
                    lower_bounds = lower, upper_bounds = upper, dimensions = len(lower), real_indices = [0],
                    integer_indices=[1], minimize = True, xopt=[0.70711, 0.5],
                    fopt= 0.750455)

def objective(x):
    x = np.array(x)
    return x[0] ** 2 + (x[1] - 2) ** 2
def constraints(x):
    x = np.array(x)
    g1 = x[1] - x[0] ** 2

    return [g1]


