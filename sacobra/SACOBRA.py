# -*- coding: utf-8 -*-
"""
Created on Fri Nov 24 10:25:12 2017

@author: r.dewinter
"""
import numpy as np

DRCL = np.array([0.3,0.05,0.001,0.0005,0.0])
DRCS = np.array([0.001,0.0])

def verboseprint(verbose, important, message):
    if verbose != 0:
        if verbose==2 or (verbose==1 and important):
            print(message)

def scaleRescale(xStart, originalL, originalU, newlower, newupper):
    return (newupper-newlower)*((xStart-originalL)/(originalU-originalL))+newlower

def rescaleWrapper(fn,lower,upper,newlower,newupper):
    oldfn = fn
    def newfn(x):
        x = scaleRescale(x,newlower,newupper,lower,upper)
        y = oldfn(x)
        return(y)
    return(newfn)
    
def adDRC(maxF, minF):
    FRange = maxF-minF
    if FRange>1e+3:
        print(FRange,'=FR is large, XI is set to short DRC')
        DRC = DRCS
    else:
        DRC = DRCL
        print(FRange,'=FR is not large, XI is set ot Long DRC')
    return DRC

def plog(y, pShift=0.0):
    if y-pShift>=0:
        ret = np.log(1+y-pShift)
    else:
        ret = -np.log(1-(y-pShift))
    return ret

def plogReverse(y,pShift=0):
    if y>0:
        ret = np.exp(y)-1+pShift
    else:
        ret = pShift+1-(1/np.exp(y))
    return ret
            
def adFit(cobra,ind):
    maxF = max(cobra['Fres'])
    minF = min(cobra['Fres'])
    FRange = maxF-minF
    
    def transfFunc(x,pShift=0.0):
        y = plog(x, pShift=pShift)
        return(y)    
    #If cobra['online'] PLOG is true then the desciosn to do the plog transfomation or not 
    #is being made in every iteration according to the p-effect otherwise the decision is made once accoridng to the FRange value 
    
    if 'onlinePlog' in cobra['sac'] and cobra['sac']['onlinePlog']:
        pShift = 0
        if cobra['pEffect'] > 1:
            Fres = np.array([transfFunc(fres, pShift) for fres in cobra['Fres'][ind]])
            if 'Plog' in cobra:
                cobra['Plog'].append(True)
            else:
                cobra['Plog'] = [True]
        else:
            Fres = cobra['Fres'][ind]
            pShift = None
            if 'Plog' in cobra:
                cobra['Plog'].append(False)
            else:
                cobra['Plog'] = [False]
    else:
        if FRange > cobra['sac']['TFRange']:
            if cobra['sac']['adaptivePLOG']:
                pShift = np.mean([cobra['fbest'],0])
            else:
                pShift = 0
            Fres = np.array([transfFunc(fres, pShift) for fres in cobra['Fres'][ind]])
            cobra['PLOG'] = [True]
        else:
            Fres = cobra['Fres'][ind]
            pShift = None
            cobra['PLOG'] = [False]
    if 'pShift' in cobra:
        cobra['pShift'].append(pShift)
    else: 
        cobra['pShift'] = [pShift]
    cobra['SurrogateInput'] = Fres
    return(cobra)


def RandomStart(cobra):
    seedn = cobra['cobraSeed']+len(cobra['Fres'])
    np.random.seed(seedn)
    anewrand = np.random.rand(1)
    
    diff = cobra['sac']['RSmax']-cobra['sac']['RSmin']
    
    if cobra['sac']['RStype']=='SIGMOID':
        tanh = np.tanh(-(len(cobra['A'])-(cobra['initDesPoints']+15)))
        randomnessTemp = (diff/2)*tanh+(diff/2)+cobra['sac']['RSmin']
    elif cobra['sac']['RStype']=='CONSTANT':
        randomnessTemp = diff/2
        
    if anewrand < randomnessTemp or cobra['progressCount'] >= cobra['sac']['Cs']:
        xStart = np.random.uniform(low=cobra['lower'], high = cobra['upper'], size=len(cobra['xbest']))
        cobra['progressCount'] = 0
    else:
        xStart = cobra['xbest']
    
    cobra['xStart'] = xStart
    return(cobra)
  
def inverseRescale(x,cobra):
    z = scaleRescale(x, cobra['lower'], cobra['upper'], cobra['originalL'], cobra['originalU'])
    return(z)
              
def getXbest(cobra):
    xbest = cobra['xbest']
    if cobra['rescale']:
        xbest = inverseRescale(xbest,cobra)
    return(xbest)
    
def getFbest(cobra):
    return(cobra['originalfn'](getXbest(cobra))[0])
        
            