# -*- coding: utf-8 -*-
"""
Created on Tue Jan 30 17:29:19 2018

@author: r.dewinter
"""
import numpy as np
from g_problem import GProblem


def problem():
    lower = np.array([0,0,-0.55,-0.55])
    upper = np.array([1200,1200,0.55,0.55])
    integer_indices=[1]
    real_indices =  np.setdiff1d(np.arange(len(lower)), integer_indices)
    return GProblem("G05", problem_type_mixint = True, objective = objective, constraints = constraints, 
                    lower_bounds = lower, upper_bounds = upper, dimensions = len(lower), real_indices = real_indices,
                    integer_indices=integer_indices, minimize = True, xopt = [679.0453, 1020.067,0.1188764, -0.3962336], fopt = 5126.4981, x_ticks = 100, x_ticks_rotate = True)

def objective(x):
    x = np.array(x)
    return 3*x[0] + 0.000001*x[0]**3 + 2*x[1] + (0.000002/3)*(x[1]**3)

def constraints(x):
    x = np.array(x)
    g1 = -x[3] + x[2] - 0.55
    g2 = -x[2] + x[3] - 0.55

    # Equality constraints
    h3 = 1000*np.sin(-x[2] - 0.25) + 1000*np.sin(-x[3] - 0.25) + 894.8 - x[0]
    h4 = 1000*np.sin(x[2] - 0.25) + 1000*np.sin(x[2] - x[3] - 0.25) + 894.8 - x[1]
    h5 = 1000*np.sin(x[3] - 0.25) + 1000*np.sin(x[3] - x[2] - 0.25) + 1294.8

    return [g1, g2, h3, h4, h5]