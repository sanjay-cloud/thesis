import numpy as np
from g_problem import GProblem

def problem():
    lower = np.zeros(13)
    upper = [1] * 9 + [100] * 3 + [1]
    integer_indices=[9,10,11]
    real_indices =  np.setdiff1d(np.arange(len(lower)), integer_indices)
    return GProblem("G01", problem_type_mixint = True, objective = objective, constraints = constraints, dimensions = len(upper),
                    lower_bounds = lower, upper_bounds = np.array(upper), 
                    integer_indices=integer_indices, real_indices = [0,1,2,3,4,5,6,7,11,12], 
                    minimize = True, fopt = -15, xopt = [1,1,1,1,1,1,1,1,1,3,3,3,1], x_ticks=10)

def objective(x):
    x = np.array(x)
    return np.sum(5*x[:4])-(5*np.sum(x[:4]**2))-(np.sum(x[4:13]))

def constraints(x):
    x = np.array(x)
    g1 = 2*x[0]+2*x[1]+x[9]+x[10] - 10
    g2 = 2*x[0]+2*x[2]+x[9]+x[11] - 10
    g3 = 2*x[1]+2*x[2]+x[10]+x[11] - 10
    
    g4 = -8*x[0]+x[9]
    g5 = -8*x[1]+x[10]
    g6 = -8*x[2]+x[11]
    
    g7 = -2*x[3]-x[4]+x[9]
    g8 = -2*x[5]-x[6]+x[10]
    g9 = -2*x[7]-x[8]+x[11]
    return [g1, g2, g3, g4, g5, g6, g7, g8, g9]

