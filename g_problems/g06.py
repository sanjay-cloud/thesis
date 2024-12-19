# -*- coding: utf-8 -*-
"""
Created on Tue Jan 30 17:29:19 2018

@author: r.dewinter
"""
import numpy as np
from g_problem import GProblem


def problem():
    lower = np.array([13,0])
    upper = np.array([100,100])
    integer_indices=[0,1]
    real_indices =  np.setdiff1d(np.arange(len(lower)), integer_indices)

    return GProblem("G06", problem_type_mixint = False, objective = objective, constraints = constraints, 
                    lower_bounds = lower, upper_bounds = upper, dimensions = len(lower),
                    integer_indices=integer_indices, real_indices = [1], minimize = True, xopt = [14.095,0.84296], fopt =  -6961.81381, x_ticks=10)

def objective(x):
    x = np.array(x)
    return (x[0] - 10)**3 + (x[1] - 20)**3

def constraints(x):
    x = np.array(x)
    g1 = -(x[0] - 5)**2 - (x[1] - 5)**2 + 100
    g2 = -82.81 + (x[0] - 6)**2 + (x[1] - 5)**2

    return [g1,g2]

