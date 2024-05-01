# -*- coding: utf-8 -*-
"""
Created on Tue Jan  9 08:57:29 2018

@author: r.dewinter
"""
import numpy as np

def BKF01(x):
    f1 = 4*(x[0]**2)+4*(x[1]**2)
    f2 = (x[0]-5)**2 + (x[1]-5)**2
    g1 = (x[0]-5)**2 + x[1]**2 - 25
    g2 = -1*((x[0]-8)**2 + (x[1]+3)**2 - 7.7)
    g3 = f1 - 140
    g4 = f2 - 50
    return np.array([f1+f2,g1,g2,g3,g4])

#minSeen = np.inf
#for x in np.linspace(0,5,50):
#    for y in np.linspace(0,3,50):
#        if BKF01([x,y])[0] < minSeen:
#            minSeen = BKF01([x,y])[0]