import numpy as np
from g_problem import GProblem

def problem():
    lower = np.array([-10,-10,-10,-10,-10,-10,-10])
    upper = np.array([10,10,10,10,10,10,10])
    integer_indices = [5,6]
    real_indices =  np.setdiff1d(np.arange(len(lower)), integer_indices)

    return GProblem("G09", problem_type_mixint = True, objective = objective, constraints = constraints, 
                    lower_bounds = lower, upper_bounds = upper, dimensions = len(lower), real_indices = real_indices,
                    integer_indices=integer_indices, minimize = True, xopt = [2.330499,1.951372,-0.4775414,4.365726,-0.6244870,1.038131,1.594227],
                    fopt = 680.6300573, x_ticks_rotate = True)

def objective(x):
    x = np.array(x)
    return (x[0] - 10)**2 + 5*(x[1] - 12)**2 + x[2]**4 + 3*(x[3] - 11)**2 + 10*x[4]**6 + 7*x[5]**2 + x[6]**4 - 4*x[5]*x[6] - 10*x[5] - 8*x[6]

def constraints(x):
    x = np.array(x)
    g1 = -127 + 2*x[0]**2 + 3*x[1]**4 + x[2] + 4*x[3]**2 + 5*x[4]
    g2 = -282 + 7*x[0] + 3*x[1] + 10*x[2]**2 + x[3] - x[4]
    g3 = -196 + 23*x[0] + x[1]**2 + 6*x[5]**2 - 8*x[6]
    g4 = 4*x[0]**2 + x[1]**2 - 3*x[0]*x[1] + 2*x[2]**2 + 5*x[5] - 11*x[6]

    return [g1, g2, g3, g4]


