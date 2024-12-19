# -*- coding: utf-8 -*-
"""
Created on Mon Jan  8 15:55:23 2018

@author: r.dewinter
"""
import numpy as np
from g_problem import GProblem


def problem():
    lower = np.array([0, 0, 0, 0, 0, 0])
    upper = np.array([10, 10, 10, 10, 10, 10])
    integer_indices=[2,5]
    real_indices =  np.setdiff1d(np.arange(len(lower)), integer_indices)

    return GProblem("G02", problem_type_mixint = True,  objective = objective, constraints = constraints, 
                    lower_bounds = lower, upper_bounds = upper, dimensions = len(lower),
                    integer_indices=integer_indices, real_indices=real_indices, minimize = True, fopt=0.03323918, xopt = [4.69338849e+00, 7.79506815e+00, 4.84204651e-05, 8.04332048e+00, 4.90557079e+00, 8.08491255e+00])


def objective(x):
    x = np.array(x)
    return np.abs((np.sum(np.cos(x)**4) - 2 * np.prod(np.cos(x)**2))/  np.sqrt(np.sum(np.arange(1,len(x)+1) * x**2)))
    

def constraints(x):
    x = np.array(x)
    g1 = 0.75 - np.prod(x)
    g2 = np.sum(x) - 7.5 * len(x)
    return [g1, g2]