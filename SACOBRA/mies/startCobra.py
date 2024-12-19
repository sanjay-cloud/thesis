# -*- coding: utf-8 -*-
"""
Created on Mon Nov 27 10:26:02 2017

@author: r.dewinter
"""
import sys
import os
sys.path.append(os.path.abspath('../../g_problems'))

import g01
import g02
import g03
import g04
import g05
import g06
import g07
import g08
import g09
import g10
import g11

# from G02_2 import G02_2
# from BKF01 import BKF01
# from G05 import G05

from cobraInit import cobraInit
from cobraPhaseII import cobraPhaseII

from SACOBRA import getXbest
from SACOBRA import getFbest

import numpy as np


# fn = G02
# fName="G02"
# nConstraints = 2
# lower = np.array([-1.5, -0.5])
# upper = np.array([1.5, 2.5])
# xStart = lower+np.random.rand(nConstraints)*upper

# fn = G02
# fName="G02_2"
# nConstraints = 2
# lower = np.array([-1.5, -0.5])
# upper = np.array([1.5, 2.5])
# xStart = lower+np.random.rand(2)*upper


#fn = BKF01
#fName="BKF01"
#nConstraints = 4
#lower = np.array([0, 0, 0, 0])
#upper = np.array([5, 3, 140, 50])
#xStart = lower+np.random.rand(nConstraints)*upper

# problems = [g01, g02, g03, g04, g05, g06, g07, g08, g09, g10, g11]


problems = [g01, g07, g09]
cobraSeeds = [0,5,10,50,100]
feval = 500
seqFeval = 1000

for problem in problems:
    problem = problem.problem()
    fn = problem.combined_obj
    fName = problem.fn()
    upper = problem.upper_bounds()
    lower = problem.lower_bounds()
    d = len(upper)
    xStart = lower+np.random.rand(len(upper))*upper    
    nConstraints = len(problem.constraints(xStart))
    integer_indices = problem.integer_indices()
    for seed in cobraSeeds:
        np.random.seed(seed)
        cobra = cobraInit(xStart, fn, fName, lower, upper, nConstraints, feval=feval, seqFeval=seqFeval, initDesPoints=3*d, DOSAC=1, cobraSeed=seed, integer_indices= integer_indices)

        cobra = cobraPhaseII(cobra)

print(getXbest(cobra))
print(getFbest(cobra))