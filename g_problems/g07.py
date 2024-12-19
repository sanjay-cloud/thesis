# -*- coding: utf-8 -*-
"""
Created on Tue Jan 30 17:29:19 2018

@author: r.dewinter
"""
import numpy as np
from g_problem import GProblem


def problem():
    lower = np.array([-10,-10,-10,-10,-10,-10,-10,-10,-10,-10])
    upper = np.array([10,10,10,10,10,10,10,10,10,10])
    integer_indices=[3,4,9]
    real_indices =  np.setdiff1d(np.arange(len(lower)), integer_indices)
    return GProblem("G07", problem_type_mixint = True, objective = objective, constraints = constraints, 
                    lower_bounds = lower, upper_bounds = upper, dimensions = len(lower), real_indices = real_indices,
                    integer_indices=integer_indices, minimize = True, xopt = [2.171906, 2.363683,8.773926, 5.095984,0.9906548,
                                                                        1.430574, 1.321644,0.828726, 8.280092 ,8.375937] ,
                                                                        fopt = 24.3062091, x_ticks_rotate = True)

def objective(x):
    x = np.array(x)
    return x[0]**2 + x[1]**2 + x[0]*x[1] - 14*x[0] - 16*x[1] + (x[2] - 10)**2 + 4*(x[3] - 5)**2 + (x[4] - 3)**2 + 2*(x[5] - 1)**2 + 5*x[6]**2 + 7*(x[7] - 11)**2 + 2*(x[8] - 10)**2 + (x[9] - 7)**2 + 45

def constraints(x):
    x = np.array(x)
    g1 = 4*x[0] + 5*x[1] - 3*x[6] + 9*x[7] - 105
    g2 = 10*x[0] - 8*x[1] - 17*x[6] + 2*x[7]
    g3 = -8*x[0] + 2*x[1] + 5*x[8] - 2*x[9] - 12
    g4 = 3*(x[0] - 2)**2 + 4*(x[1] - 3)**2 + 2*x[2]**2 - 7*x[3] - 120
    g5 = 5*x[0]**2 + 8*x[1] + (x[2] - 6)**2 - 2*x[3] - 40
    g6 = 0.5*(x[0] - 8)**2 + 2*(x[1] - 4)**2 + 3*x[4]**2 - x[5] - 30
    g7 = x[0]**2 + 2*(x[1] - 2)**2 - 2*x[0]*x[1] + 14*x[4] - 6*x[5]

    return [g1, g2, g3, g4, g5, g6, g7]

