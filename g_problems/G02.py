# -*- coding: utf-8 -*-
"""
Created on Mon Jan  8 15:55:23 2018

@author: r.dewinter
"""
import numpy as np

def G02(x):
    obj = (1-x[0])**2 + 100*(x[1]-x[0]**2)**2
    g1 = (x[0]-1)**3-x[1]+1
    g2 = x[0]+x[1]-2
    res = np.array([obj,g1,g2])
    return(res)