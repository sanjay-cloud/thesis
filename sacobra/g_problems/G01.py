# -*- coding: utf-8 -*-
"""
Created on Fri Nov 24 09:31:35 2017

@author: r.dewinter
"""
import numpy as np

def G01(x):
    obj = np.sum(5*x[:4])-(5*np.sum(x[:4]**2))-(np.sum(x[4:13]))
    g1 = 2*x[0]+2*x[1]+x[9]+x[10] - 10
    g2 = 2*x[0]+2*x[2]+x[9]+x[11] - 10
    g3 = 2*x[1]+2*x[2]+x[10]+x[11] - 10
    
    g4 = -8*x[0]+x[9]
    g5 = -8*x[1]+x[10]
    g6 = -8*x[2]+x[11]
    
    g7 = -2*x[3]-x[4]+x[9]
    g8 = -2*x[5]-x[6]+x[10]
    g9 = -2*x[7]-x[8]+x[11]
    
    return(np.array([obj,g1,g2,g3,g4,g5,g6,g7,g8,g9]))