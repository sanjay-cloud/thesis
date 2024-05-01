# -*- coding: utf-8 -*-
"""
Created on Mon Nov 11 15:58:15 2019

@author: r.dewinter
"""
import numpy as np


def cobraPhaseI(cobra):
    print("No feasible point available in initial sample")
    print("Phase I started")
    
    phase = "PHASE1"
    fn = cobra['fn']
    dimension = cobra['dimension']
    nConstraints = cobra['nConstraints']
    EPS = np.array([cobra['epsilonInit']]*cobra['nConstraints'])
    n = len(cobra['A'])
    iteration = cobra['iteration']
    feasibleSolutionExists = 0 in cobra['numViol']
    predY = np.empty(cobra['initDesPoints'])
    predY[:] = np.nan
    cobra['predC'] = np.empty((cobra['initDesPoints'], cobra['nConstraints'])) # matrix to store predicted constraint values
    cobra['predC'][:] = np.nan
    constraintPrediction = None
    optimizerConvergence = np.ones(cobra['initDesPoints'])
    # maybe finish this, but it does not seem logic since default poarameters tell that this will never be used it is not used.    
    