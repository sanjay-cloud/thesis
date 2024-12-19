# -*- coding: utf-8 -*-
"""
Created on Mon Jan  8 15:55:23 2018

@author: r.dewinter
"""
import numpy as np
from g_problem import GProblem


def problem():
    upper = np.array([102,45,45,45,45])
    lower = np.array([78,33,27,27,27])
    integer_indices=[0,1]
    real_indices =  np.setdiff1d(np.arange(len(lower)), integer_indices)

    return GProblem("G04", problem_type_mixint = True, objective = objective, constraints = constraints, 
                    lower_bounds = lower, upper_bounds = upper, dimensions=len(lower),
                    integer_indices=integer_indices, real_indices = real_indices, minimize = True, fopt= -30665.5, 
                    xopt=[78.0,33.0,29.995,45.0,36.776], x_ticks = 2, x_ticks_rotate = True)


def objective(x):
    x = np.array(x)
    return 5.3578547 * x[2]**2 + 0.8356891 * x[0] * x[4] + 37.293239 * x[0] - 40792.141

def constraints(x):
    x = np.array(x)
    g1 = 85.334407 + 0.0056858 * x[1] * x[4] + 0.0006262 * x[0] * x[3] - 0.0022053 * x[2] * x[4] - 92
    g2 = -85.334407 - 0.0056858 * x[1] * x[4] - 0.0006262 * x[0] * x[3] + 0.0022053 * x[2] * x[4]
    g3 = 80.51249 + 0.0071317 * x[1] * x[4] + 0.0029955 * x[0] * x[1] + 0.0021813 * x[2]**2 - 110
    g4 = -80.51249 - 0.0071317 * x[1] * x[4] - 0.0029955 * x[0] * x[1] - 0.0021813 * x[2]**2 + 90
    g5 = 9.300961 + 0.0047026 * x[2] * x[4] + 0.0012547 * x[0] * x[2] + 0.0019085 * x[2] * x[3] - 25
    g6 = -9.300961 - 0.0047026 * x[2] * x[4] - 0.0012547 * x[0] * x[2] - 0.0019085 * x[2] * x[3] + 20
    
    return [g1, g2, g3, g4, g5, g6]