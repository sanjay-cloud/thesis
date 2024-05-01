# -*- coding: utf-8 -*-
"""
Created on Tue Nov 21 16:46:10 2017

@author: r.dewinter
"""
from SACOBRA import DRCL
from SACOBRA import scaleRescale
from SACOBRA import rescaleWrapper
from SACOBRA import adDRC
from SACOBRA import verboseprint
from defaultRI import defaultRI
from defaultTR import defaultTR
from defaultSAC import defaultSAC
from defaultSAC import setOpts
from lhs import lhs 
from transformLHS import transformLHS


import warnings
import numpy as np

def cobraInit(xStart, fn, fName, lower, upper, nConstraints, feval, 
              initDesign="LHS", 
              initDesPoints=None, 
              initDesOptP=None, 
              initBias=0.005,
              seqOptimizer="COBYLA", 
              seqFeval=1000, 
              seqTol=1e-6, 
              penaF=[3.0, 1.7, 3e5], 
              squaresF=True, 
              squaresC=True, 
              conTol=0.0,
              constraintHandling="DEFAULT",
              sigmaD=[3.0,2.0,100], 
              repairInfeas=False, 
              ri=defaultRI(),
              DOSAC=1, 
              sac=None, 
              epsilonInit=None, 
              epsilonMax=None, 
              solu=None,
              saveIntermediate=False, 
              saveSurrogates=False, 
              RBFmodel="cubic", 
              RBFwidth=-1,
              GaussRule="One",
              RBFrho=0.0,
              skipPhaseI=True,
              trueFuncForSurrogates=False, 
              drFactor=1.0, 
              XI=DRCL,
              rescale=True,
              newlower=-1,
              newupper=1,  
              DEBUG_XI=False, 
              SKIP_HIGH=False, 
              DEBUG_RBF=False,
              TrustRegion=False,
              TRlist=defaultTR(),
              verbose=1,
              verboseIter=10,
              cobraSeed=1):
    
    if initDesPoints is None:
        initDesPoints = 2*len(xStart)+1
    if sac is None:
        sac = defaultSAC(DOSAC)
        
    originalfn = fn
    originalL = lower
    originalU = upper
    phase = 'init'
    
    dimension = len(xStart) #number of parameters
    
    if rescale:
        lower = np.array([newlower]*dimension)
        upper = np.array([newupper]*dimension)
        xStart = scaleRescale(xStart, originalL, originalU, newlower, newupper)
        fn = rescaleWrapper(fn,originalL,originalU,newlower,newupper)
    l = min(upper-lower)
    
    if epsilonInit is None:
        epsilonInit = 0.005*l
    if epsilonMax is None:
        epsilonMax = 2*0.005*l
    if initDesOptP is None:
        initDesOptP = initDesPoints
        
    if initDesPoints>=feval:
        raise ValueError('feval < initDesPoints is not true')
    #if len(fn(xStart))!=(nConstraints+1): #DOM!!!!!!!! BECAUSE OF NUMBER OF FN
    #    raise ValueError('Numer of constraints and objectives computed with fn is not equal to number of constraints+1')
    
    np.random.seed(cobraSeed)
    dimension = len(xStart)
    iteration = 0
    
    I = np.empty((1,1))
    I[:] = np.NaN
    Gres = np.empty((1,1))
    Gres[:] = np.NaN
    Fres = []
    
    if initDesign == 'RANDOM':
        np.random.seed(cobraSeed)
        I = np.random.uniform(low=lower,high=upper,size=(initDesPoints-1,dimension))
        I = np.vstack((I,xStart))
        I = np.clip(I,lower,upper)
        
        randomResults = randomResultsFactory(I,fn,nConstraints)
        Gres = randomResults[1:].T
        Fres = randomResults[0]
    elif initDesign =='LHS':
        np.random.seed(cobraSeed)
        I = lhs(dimension, samples=initDesPoints-1, criterion="center", iterations=5)
        I = transformLHS(I, lower, upper)
        I = np.vstack((I,xStart))
        I = np.clip(I,lower,upper)

        randomResults = randomResultsFactory(I,fn,nConstraints)
        Gres = randomResults[1:].T
        Fres = randomResults[0]        
    else:
        raise ValueError('not yet implemented or invalid init design')
    
    Tfeas = np.floor(2*np.sqrt(dimension)) # The threshhold parameter for the number of consecutive iterations that yield feasible solution before the margin is reduced
    Tinfeas = np.floor(2*np.sqrt(dimension)) # The threshold parameter for the number of consecutive iterations that yield infeasible solutions before the margin is increased
    
    fe = len(I) #number of function evaluations
    fbest = []
    xbest = []
    fbestArray = []
    xbestArray = np.empty((1,1))
    xbestArray[:] = np.NaN
    
    numViol = np.sum(Gres>0,axis=1)
    
    maxViol = np.max([np.zeros(len(Gres)),np.max(Gres, axis=1)],axis=0)
    
    A = I #contains all evaluated points
    n = len(A)
    
    # determining best feasible objective value (fbest) so far and best point (xbest)
    if 0 in numViol:
        fbestI = np.argmin(Fres[numViol==0])
        fbest = Fres[fbestI]
        xbest = I[fbestI]
        ibest = fbestI
    else:
        minNumIndex = np.where(numViol==min(numViol))[0]
        FresMin = Fres[minNumIndex]
        ind = np.argmin(FresMin)
        index = minNumIndex[ind]
        
        fbest = Fres[index]
        xbest = A[index]
        ibest = index
    
    fbestArray = np.array([fbest]*initDesPoints)
    xbestArray = np.tile(xbest,(n,1))
    
    print('start run with seed',cobraSeed)
    print('Initialization is done')
    
    cobra = dict()
    cobra['CONSTRAINED'] = True
    cobra['fn'] = fn
    cobra['xStart'] = xbest
    cobra['fName'] = fName
    cobra['dimension'] = dimension
    cobra['nConstraints'] = nConstraints
    cobra['lower'] = lower
    cobra['upper'] = upper
    cobra['newlower'] = newlower
    cobra['newupper'] = newupper
    cobra['originalL'] = originalL
    cobra['originalU'] = originalU
    cobra['originalfn'] = originalfn
    cobra['rescale'] = rescale
    cobra['feval'] = feval
    cobra['A'] = A
    cobra['fbestArray'] = fbestArray
    cobra['xbestArray'] = xbestArray
    cobra['xbest'] = xbest
    cobra['fbest'] = fbest
    cobra['ibest'] = ibest
    cobra['Fres'] = Fres
    cobra['Gres'] = Gres
    cobra['numViol'] = numViol
    cobra['maxViol'] = maxViol
    cobra['epsilonInit'] = epsilonInit
    cobra['epsilonMax'] = epsilonMax
    cobra['XI'] = XI
    cobra['drFactor'] = drFactor
    cobra['Tinfeas'] = Tinfeas
    cobra['Tfeas'] = Tfeas
    cobra['iteration'] = iteration
    cobra['initDesPoints'] = initDesPoints
    cobra['seqOptimizer'] = seqOptimizer
    cobra['ptail'] = True
    cobra['seqFeval'] = seqFeval
    cobra['seqTol'] = seqTol
    cobra['penaF'] = penaF
    cobra['sigmaD'] = sigmaD
    cobra['squaresF'] = squaresF
    cobra['squaresC'] = squaresC
    cobra['cobraSeed'] = cobraSeed
    cobra['conTol'] = conTol
    cobra['constraintHandling'] = constraintHandling
    cobra['l'] = l
    cobra['repairInfeas'] = repairInfeas
    cobra['ri'] = ri
    cobra['fe'] = fe
    cobra['saveIntermediate'] = saveIntermediate
    cobra['saveSurrogates'] = saveSurrogates
    cobra['RBFmodel'] = RBFmodel
    cobra['RBFwidth'] = RBFwidth
    cobra['RBFrho'] = RBFrho
    cobra['RULE'] = GaussRule
    cobra['widthFactor'] = 1.0
    cobra['skipPhaseI'] = skipPhaseI
    cobra['trueFuncForSurrogates'] = trueFuncForSurrogates
    cobra['solu'] = solu
    cobra['TrustRegion'] = TrustRegion
    cobra['TRlist'] = TRlist
    cobra['radi'] = np.array([TRlist['radiInit']]*initDesPoints)
    cobra['sac'] = sac
    cobra['TRDONE'] = [False]*initDesPoints
    cobra['TRind'] = []
    cobra['refinedX'] = []
    cobra['fCount'] = 0
    cobra['sCount'] = 0
    cobra['DOSAC'] = DOSAC
    cobra['PLOG'] = False
    cobra['pShift'] = [0]
    cobra['pEffect'] = None
    cobra['progressCount'] = 0
    cobra['DEBUG_XI'] = DEBUG_XI
    cobra['DEBUG_RBF'] = DEBUG_RBF
    cobra['SKIP_HIGH'] = SKIP_HIGH
    cobra['verbose'] = verbose
    cobra['verboseIter'] = verboseIter
    cobra['phase'] = [phase]*initDesPoints
        
    cobra['ri'] = setOpts(cobra['ri'],defaultRI())
    
    if DOSAC>0:
        verboseprint(cobra['verbose'], False, 'Parameter and function adjustment phase')
        cobra['sac'] = setOpts(cobra['sac'], defaultSAC(DOSAC))
        if 'sac' in cobra:
            cobra['pEffect'] = cobra['sac']['pEffectInit']
        else:
            cobra['pEffect'] = 0
        cobra['TRlist'] = setOpts(cobra['TRlist'], defaultTR())
        
        if 'sac' in cobra:
            if 'aDRC' in cobra['sac']:
                if len(XI)!=len(DRCL):
                    warnings.warn('Xi is differnt from default DRCL but sac[aDRC] is true, so XI will be set by automatic DRC adjustment' ,DeprecationWarning)
                elif np.any(XI!=DRCL):
                    warnings.warn('Xi is differnt from default DRCL but sac[aDRC] is true, so XI will be set by automatic DRC adjustment' ,DeprecationWarning)
            verboseprint(cobra['verbose'], False, "adjusting DRC")
            DRC = adDRC(np.max(cobra['Fres']),np.min(cobra['Fres']))
            cobra['XI'] = DRC
    return(cobra)

def randomResultsFactory(I,fn,nConstraints):
    result = np.empty((nConstraints+1,len(I)))
    for i in range(len(I)):
        result[:,i] = fn(I[i])
    return result
        


            
