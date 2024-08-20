# -*- coding: utf-8 -*-
"""
Created on Mon Nov 27 10:26:02 2017

@author: r.dewinter
"""
from G01 import G01
from G02 import G02
from G02_2 import G02_2
from BKF01 import BKF01
from G05 import G05

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
 
cobraSeed = 1
np.random.seed(cobraSeed)

fn = G01
fName = 'G01'
upper = np.array([1,1,1,1,1,1,1,1,1,100,100,100,1])
lower = np.zeros(len(upper))
nConstraints = 9
np.random.seed(1)
d = 13
xStart = lower+np.random.rand(len(upper))*upper
#xStart = np.array([0.2655087,0.3721239,0.5728534,0.9082078,0.2016819,0.8983897,0.9446753,0.6607978,0.6291140,6.1786270,20.5974575,17.6556753,0.6870228])

# fn = G05
# fName= 'G05'
# lower = np.array([0,0,-0.55,-0.55])
# upper = np.array([1200,1200,0.55,0.55])
# nConstraints=8
# d=4
# xStart = lower+np.random.rand(len(upper))*upper
feval = 150
seqFeval = 1000
cobra = cobraInit(xStart, fn, fName, lower, upper, nConstraints, feval=feval, seqFeval=seqFeval, initDesPoints=3*d, DOSAC=1, cobraSeed=1)

cobra = cobraPhaseII(cobra)

print(getXbest(cobra))
print(getFbest(cobra))