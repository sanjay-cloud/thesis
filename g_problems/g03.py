# -*- coding: utf-8 -*-
"""
Created on Mon Jan  8 15:55:23 2018

@author: r.dewinter
"""
import numpy as np
from g_problem import GProblem


def problem():
    upper = np.array([1,1,1,1,1])
    lower = np.zeros(len(upper))
    integer_indices=[0,3]
    real_indices =  np.setdiff1d(np.arange(len(lower)), integer_indices)

    return GProblem("G03", problem_type_mixint = True, objective = objective, constraints = constraints, 
                    lower_bounds = lower, upper_bounds = upper, dimensions = len(upper),
                    integer_indices=integer_indices, real_indices=real_indices, minimize = False, fopt = 0.990866, 
                    xopt = [0.52745941, 0.5 ,  0.33279939 ,0.53295108, 0.5  ])


def objective(x):
    x = np.array(x)
    n = len(x)

    return (np.sqrt(n) ** n) * np.prod(x)

def constraints(x):
    x = np.array(x)
    g1 = np.sum(np.square(x)) - 1
    return [g1]