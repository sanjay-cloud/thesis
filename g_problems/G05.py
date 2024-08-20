# -*- coding: utf-8 -*-
"""
Created on Tue Jan 30 17:29:19 2018

@author: r.dewinter
"""
import numpy as np

def G05(x):
    print("yes")
    obj =  3*x[0] + 0.000001*x[0]**3 + 2*x[1] + (0.000002/3)*x[1]**3
    g1 = -x[3] + x[2] - 0.55
    g2 = -x[2] + x[3] - 0.55
    
    h31 = -1*(1000*np.sin(-x[2]-0.25) + 1000*np.sin(-x[3]-0.25) + 894.8 - x[0] + 0.0001)
    h32 = 1000*np.sin(-x[2]-0.25) + 1000*np.sin(-x[3]-0.25) + 894.8 - x[0] - 0.0001
    h41 = -1*(1000*np.sin(x[2]-0.25) + 1000*np.sin(x[2]-x[3]-0.25) + 894.8 - x[1] + 0.0001)
    h42 = 1000*np.sin(x[2]-0.25) + 1000*np.sin(x[2]-x[3]-0.25) + 894.8 - x[1] - 0.0001
    h51 = -1*(1000*np.sin(x[3]-0.25) + 1000*np.sin(x[3]-x[2]-0.25) + 1294.8 + 0.0001)
    h52 = 1000*np.sin(x[3]-0.25) + 1000*np.sin(x[3]-x[2]-0.25) + 1294.8 - 0.0001

    return(np.array([obj,g1,g2,h31,h32,h41,h42,h51,h52]))