# -*- coding: utf-8 -*-
"""
Created on Wed Jan 10 14:52:00 2018

@author: r.dewinter
"""

import numpy as np

def G02_2(x):
    obj = (1-x[0])**2 + 100*(x[1]-x[0]**2)**2
    g1 = (x[0]-1)**3-x[1]+1
    g2 = x[0]+x[1]-2
    g3 = x[0]-1.00000001
    g4 = -x[0]+0.99999999
    res = np.array([obj,g1,g2,g3,g4])
    return(res)