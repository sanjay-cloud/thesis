# -*- coding: utf-8 -*-
"""
Created on Mon Nov 27 10:48:42 2017

@author: r.dewinter
"""

from SACOBRA import verboseprint
from SACOBRA import adFit
from SACOBRA import plog
from SACOBRA import RandomStart
from SACOBRA import plogReverse
from SACOBRA import getXbest

from RbfInter import trainCubicRBF
from RbfInter import predictRBFinter
from RbfInter import distLine
from RbfInter import interpRBF


from scipy import optimize
import numpy as np
import time
import warnings
import copy
import json

##########################################################################################################
# Some rules about the COBRA-II-code:
#
# - cobra$df contains one row for each iteration, including initial points
#   cobra$df2 contains one row for each phase-II iteration only
# - cobra$PLOG is set by adFit (SACOBRA.R) and adFit is called at the start of each iteration
#   (see trainSurrogates here in phase II)
# - cobra$solu is always in original input space. But when it is used (for diagnostics) in 
#   updateSaveCobra, then a local copy \code{solu} is made, and - if cobra$rescale==TRUE - 
#   \code{solu} is transformed to the rescaled space.
# - the rescaled space has the bounds [rep(cobra$newlower,d),rep(cobra$newupper,d)], (usually -1 and 1) 
# - cobra$xbest,cobra$fbest,cobra$ibest refer always to the same infill point (same iteration).
# - What is cobra$fbest before a feasible infill is found? - See \code{updateSaveCobra}:
#   If the new infill point is feasible, take its fn[1]-value (of course!). If it is not feasible,
#   leave it at the setting from cobraInitial.R#400: from all points of the initial design
#   with minimum number of violated constraints, take the one with smallest Fres.
# 
##########################################################################################################
def cobraPhaseII(cobra):
    '''
    Improve the feasible solution by searching new infill points
    Improve the feasible solution using the COBRA optimizer phase II
    by searching new infill points with the help of RBF surrogate models. 
    May be even called if no feasible solution is found yet, then phase II will try to find
    feasible solutions.
    @param cobra an object of class COBRA, this is a (long) list containing all settings
        from cobraInit
    @return cobra, an object of class COBRA from cobraInit, 
      enhanced here by the following elements (among others):
        fn function returning an (m+1)-vector c(objective,g1,...,gm). This
              function may be a rescaled and plog-transformed version of the original fn 
              passed into cobraInit. The original fn is in 
              cobra$originalFn. 
        df  data frame with summary of the optimization run (see below)
        df2  data frame with additional summary information (see below)
        A (feval x dim)-matrix containing all evaluated points 
              in input space. If rescale==TRUE, all points are in rescaled input space. 
        Fres a vector of the objective values of all evaluated points 
        Gres a matrix of the constraint values of all evaluated points 
        fbest the best feasible objective value found 
        xbest the point in input space yielding the best feasible objective value 
        ibest the corresponding iteration number (row of cobra$df, of cobra$A)
        PLOG If TRUE, then the objective surrogate model is trained on the 
              plog-transformed objective function. 
        
     Note that cobra$Fres, cobra$fbest, cobra$fbestArray and similar contain 
     always the objective values of the orignial function cobra$fn[1]. (The surrogate models 
     may be trained on a plog-transformed version of this function.)
     
     The data frame cobra$df contains one row per iteration with columns 
        iter  
        y   true objective value Fres 
        predY  surrogate objective value. Note: The surrogate may be trained on  
              plog-transformed training data, but predY is transformed back to the original 
              objective range. NA for the initial design points.
        predSolu  surrogate objective value at best-known solution cobra$solu, if given. 
              If cobra$solu is NULL, take the current point instead. Note: The surrogate may be trained on  
              plog-transformed training data, but predSolu is transformed back to the original 
              objective range. NA for the initial design points.
        feasible  
        feasPred  
        nViolations  
        maxViolation  
        FEval  number of function evaluations in sequential optimizer. NA if it was a repair step 
        Best  ever-best feasible objective value fbest. As long as there is 
              no feasible point, take among those with minimum number of violated constraints the
              one with minimum Fres. 
        optimizer e.g. "COBYLA"  
        optimizationTime  in sec
        conv  
        seed  
     
     The data frame cobra$df2 contains one row per phase-II-iteration with columns 
        iter  
        predY  surrogate objective value. Note: The surrogate may be trained on  
              plog-transformed training data, but predY is transformed back to the original 
              objective range. NA for the initial design points.
        predVal   surrogate objective value + penalty 
        predSolu   surrogate objective value at true solution (see cobra$df$predSolu) 
        predSoluPenal   surrogate objective value + penalty at true solution (only diagnostics)
        sigmaD  
        penaF  
        XI  the DRC element used in the current iteration 
        EPS  
     
    seealso   cobraPhaseI, cobraInit
    '''
    verboseprint(cobra['verbose'],False,"There is at least one feasible point in the population or PhaseI is skipped")
    verboseprint(2, True, "PHASE II started")
    phase = 'PHASE II'
    if cobra['fbest'] is None:
        raise ValueError("cobraPhaseII: cobra['fbest'] is None!")
    if cobra['ibest'] is None:
        raise ValueError("cobraPhaseII: cobra['ibest'] is None!")
    ###########################################################################
    # STEP5:                                                                  #
    # Initializing the parameters and Initialize the margin  and counters     #
    ###########################################################################
    newErr1 = 0
    newErr2 = 0
    err1 = np.array([])
    err2 = np.array([])
    fn = cobra['fn']
    CHECKDIST = True
    Tfeas = cobra['Tfeas']
    Tinfeas = cobra['Tinfeas']
    Cfeas = 0 # Starting Counters
    Cinfeas = 0
    EPS = np.array([cobra['epsilonInit']]*cobra['nConstraints']) # Initializing margin for all constraints
    n = len(cobra['A'])
    nRepair = 0
    if n==cobra['initDesPoints']:
        predY = np.empty(cobra['initDesPoints'])
        predY[:] = np.nan # structure to store surrogate optimization results
        predVal = np.empty(cobra['initDesPoints'])
        predVal[:] = np.nan
        if cobra['nConstraints']!=0:
            cobra['predC'] = np.empty((cobra['initDesPoints'], cobra['nConstraints'])) # matrix to store predicted constraint values
            cobra['predC'][:] = np.nan
            feas = np.array([not any(cobra['Gres'][i]>0) for i in range(len(cobra['Gres']))])# feasibility of initial design
        feasPred = np.empty(cobra['initDesPoints'])
        feasPred[:] = np.nan
        optimizerConvergence = np.ones(cobra['initDesPoints']) # vector to store optimizer convergence
        cobra['optimizationTime'] = np.zeros(cobra['initDesPoints'])
        feval = np.empty(cobra['initDesPoints'])
        feval[:] = np.nan
        fbestArray = np.empty(cobra['initDesPoints'])
        fbestArray[:] = cobra['fbest']
    else:
        raise ValueError('to be implemented')
    
    constraintSurrogates = None
    fitnessSurrogate = None
    fitnessSurrogate1 = None
    fitnessSurrogate2 = None
    
    constraintPrediction = None
    penaF = cobra['penaF']
    sigmaD = cobra['sigmaD']
    cobra['important'] = False
    
    if n >= cobra['feval']:
        print("n is after Phase I equal or larger than cobra['feval']")
        
    def fitFuncPenalRBF(x):
        nonlocal fitnessSurrogate
        nonlocal ro
        def distRequirement(x, fitnessSurrogate,ro):
            ed = ro - distLine(x,fitnessSurrogate['xp'])
            violatedDist = np.where(ed>0)[0]
            return (sum(ed[violatedDist]))
        
        if np.any(np.isnan(x)):
            warnings.warn('fitFuncPenalRBF: x value is Nan, returning Inf')
            return(np.inf)
        
        y = interpRBF(x,fitnessSurrogate)
        nonlocal cobra
        nonlocal fn
        if cobra['trueFuncForSurrogates']:
            y = fn(x)[0]
        nonlocal constraintPrediction
        nonlocal constraintSurrogates
        cstr = interpRBF(x,constraintSurrogates)+EPS**2 
        constraintPrediction =cstr
        if cobra['trueFuncForSurrogates']:
            cstr = fn(x)[1:]+EPS**2
            constraintPrediction = cstr
        violatedConsraints = np.where(cstr>0)[0]
        penalty = np.sum(cstr[violatedConsraints])
        
        nonlocal sigmaD
        penalty = penalty + distRequirement(x,fitnessSurrogate,ro)*sigmaD[0]
        return(y+penalty*penaF[0])
        
    def checkDistanceReq(subMin, constraintSurrogates, sigmaD, CHECKDIST):
        if CHECKDIST:
            nonlocal ro
#            global ro
            ed = ro-distLine(subMin['x'], constraintSurrogates['xp'])
            violatedDist = np.where(ed>0)[0]
            if len(violatedDist)>0:
                if sigmaD[0]*sigmaD[1]<sigmaD[2]:
                    sigmaD[0] = sigmaD[0]*sigmaD[1]
                    verboseprint(cobra['verbose'], False, "increasing sigmaD to: "+str(sigmaD[0])+" at iteration "+str(n))
        return(sigmaD)
        
    def updateInfoAndCounters(cobra, xNew, xNewEval, newNumViol, newMaxViol, phase):
        cobra['A'] = np.vstack((cobra['A'], xNew))
        cobra['Fres'] = np.append(cobra['Fres'], xNewEval[0])
        cobra['Gres'] = np.vstack((cobra['Gres'], xNewEval[1:]))
        cobra['numViol'] = np.append(cobra['numViol'], newNumViol)
        cobra['maxViol'] = np.append(cobra['maxViol'], newMaxViol)
        cobra['phase'] = np.append(cobra['phase'], phase)
        if len(cobra['A'])%cobra['verboseIter']==0:
            cobra['important']=True
        else:
            cobra['important']=False        
        xNewIndex = len(cobra['numViol'])-1

        realXbest = getXbest(cobra)
        verboseprint(cobra['verbose'], cobra['important'], 'Best Result '+phase+' '+str(len(cobra['A']))+': objValue = '+str(cobra['fbest'])+', x='+str(realXbest)+', viol='+str(cobra['maxViol'][cobra['ibest']]))
        
        
        nonlocal Cfeas
        nonlocal Cinfeas
        if cobra['numViol'][xNewIndex]==0:
            Cfeas += 1
            Cinfeas = 0
        else:
            Cinfeas += 1
            Cfeas = 0
        return(cobra)
        
    def adjustMargins(Cfeas,Tfeas,Cinfeas,Tinfeas,EPS,epsMax):
        if Cfeas>=Tfeas:
            EPS = EPS/2
            
            verboseprint(cobra['verbose'], False, 'reducing epsilon to: '+str(EPS[0]))
            Cfeas = 0
        
        if Cinfeas >= Tinfeas:
            EPS = np.minimum(2*EPS,epsMax)
            verboseprint(cobra['verbose'], False, 'increasing epsilon to: '+str(EPS[0]))
            Cfeas = 0
        return(Cfeas, EPS)
    
    def updateSaveCobra(cobra, xNew, feas, feasPred, feval, optimizerConvergence,
                        predY, predVal, subMin, sigmaD, penaF, gama, EPS):
        xNewIndex = len(cobra['numViol'])-1
        
        if cobra['numViol'][cobra['ibest']]==0:
            if cobra['numViol'][xNewIndex]==0 and cobra['Fres'][xNewIndex]<cobra['fbest']:
                cobra['xbest'] = xNew
                cobra['fbest'] = cobra['Fres'][xNewIndex]
                cobra['ibest'] = xNewIndex
                cobra['progressCount'] = 0
        else:
            if cobra['numViol'][xNewIndex]==0:
                cobra['xbest'] = xNew
                cobra['fbest'] = cobra['Fres'][xNewIndex]
                cobra['ibest'] = xNewIndex
                cobra['progressCount'] = 0
                
        cobra['fbestArray'] = np.append(cobra['fbestArray'],cobra['fbest'])
        cobra['xbestArray'] = np.vstack((cobra['xbestArray'],cobra['xbest']))
        
        solu = cobra['solu']
        if solu is None:
            solu=subMin['x']
        else:
            raise ValueError('to be implemented')

        if cobra['saveSurrogates']:
            cobra['constraintSurrogates'] = constraintSurrogates
            cobra['fitnessSurrogate'] = fitnessSurrogate
        
        if cobra['saveIntermediate']:
            cobraResult = copy.deepcopy(cobra)
            for d in cobraResult:
                d['corr'] = 'corr'
                d['regr'] = 'regr'
                for key in d:
                    if type(d[key]) is np.ndarray:
                        d[key] = d[key].tolist()
            
            with open('gauss_model.json', 'w') as fOut:
                json.dump(cobraResult, fOut)
        
        return(cobra)
        
        
        
    def trainSurrogates(cobra):
        nonlocal printP
        verboseprint(cobra['verbose'], False, "training "+cobra['RBFmodel']+" surrogates...")

        if cobra['SKIP_HIGH']:
            quFres = np.percentile(cobra['Fres'], 90)
            ind = np.where(cobra['Fres']<=quFres)[0]
        else:
            ind = np.arange(0,len(cobra['Fres']))
        
        A = cobra['A'][ind]
        
        if cobra['sac']['aFF']:
            cobra = adFit(cobra,ind)
            Fres = cobra['SurrogateInput']
        else:
            Fres = cobra['Fres'][ind]
        
        if cobra['DOSAC']>0:
            if cobra['PLOG'][len(cobra['PLOG'])-1] and printP:
                verboseprint(cobra['verbose'], True, "PLOG transformation is done")
        printP = False
        
        Gres = cobra['Gres'][ind]
        
        if cobra['RBFmodel'] == 'cubic':
            nonlocal constraintSurrogates
            constraintSurrogates = []
            for g in Gres.T:
                gresi = np.array(g)
                constraintSurrogates.append(trainCubicRBF(A,gresi,ptail=cobra['ptail'],squares=cobra['squaresC'],rho=cobra['RBFrho']))
                
            nonlocal fitnessSurrogate
            fitnessSurrogate = trainCubicRBF(A,Fres,ptail=cobra['ptail'],squares=cobra['squaresF'],rho=cobra['RBFrho'])
        else:
            raise ValueError('To be implemented or invalid RBFmodel')
        
        if cobra['DEBUG_RBF']:
            raise ValueError('To be implemented')
            
        #possibilty to measure p-effect after every 10 iterations
        if (cobra['sac']['onlinePLOG'] and len(cobra['A'])%cobra['sac']['onlineFreqPLOG']==0) or len(cobra['A'])==cobra['initDesPoints']:
            Fres1 = cobra['Fres']
            Fres2 = np.array([plog(fre) for fre in cobra['Fres']])
            
            if cobra['RBFmodel']=='cubic':
                nonlocal fitnessSurrogate1
                fitnessSurrogate1 = trainCubicRBF(A,Fres1,ptail=cobra['ptail'],squares=cobra['squaresF'],rho=cobra['RBFrho'])
                nonlocal fitnessSurrogate2
                fitnessSurrogate2 = trainCubicRBF(A,Fres2,ptail=cobra['ptail'],squares=cobra['squaresF'],rho=cobra['RBFrho'])
            else:
                raise ValueError('To be implemented or invalid RBFmodel')
        
        return(cobra)
        
    def getPredY(xNew, fitnessSurrogate, cobra):
        predy = interpRBF(np.asmatrix(xNew), fitnessSurrogate)
        if cobra['PLOG'][len(cobra['PLOG'])-1]:
            predy = plogReverse(predy,cobra['pShift'][len(cobra['pShift'])-1])
        return(predy)
        
    def getPredY1(xNew, fitnessSurrogate, cobra):
        predy = interpRBF(np.asmatrix(xNew), fitnessSurrogate)
        return(predy)
        
    def subProb2(x):
        nonlocal cobra
        nonlocal fitnessSurrogate
#        global cobra
#        global fitnessSurrogate
        if np.any(np.isnan(x)):
            return np.inf
        if cobra['trueFuncForSurrogates']:
            y = cobra['fn'](x)[0]
        else:
            y = predictRBFinter(fitnessSurrogate,np.asmatrix(x))
        return(y[0])
    
    def getConstraintPrediction(x, constraintSurrogates, EPS=None):        
        constraintPrediction = []
        consi = 0
        for constraintSurrogate in constraintSurrogates:
            if EPS is None:
                constraintPrediction.append(interpRBF(np.asmatrix(x),constraintSurrogate))
            else:
                constraintPrediction.append(interpRBF(np.asmatrix(x),constraintSurrogate)+EPS[consi]**2)
            consi += 1
        constraintPrediction = np.array(constraintPrediction)
        return constraintPrediction
    
    def gCOBRA(x):
        nonlocal cobra
        nonlocal ro
        nonlocal constraintPrediction
        nonlocal EPS
        nonlocal constraintSurrogates
#        global cobra
#        global ro
#        global constraintPrediction
#        global EPS
#        global constraintSurrogates
        x = np.matrix(x)
        distance = distLine(x, cobra['A'])
        subC = np.maximum(ro-distance,np.zeros(len(distance)))
        h = np.sum(subC)*cobra['drFactor']

        constraintPrediction = getConstraintPrediction(x, constraintSurrogates, EPS)
        
        if np.any(np.isnan(x)):
            warnings.warn('gCOBRA: x value is NaN, returning Inf',DeprecationWarning)
            return(np.full(x.shape, np.inf))
        
        h = np.append(np.array([-1*h]), -1*constraintPrediction)   
        return(h) 

    
    ###########################################################################
    # STEP6:                                                                  #
    #  Improve the feasible point                                             #
    ###########################################################################
    
    printP = True
    
    while n < cobra['feval']:
        ##########################################################
        # STEP6.1: UPDATE RBF MODEL for fitness and constratnit  #
        ##########################################################

        #checken of dit goed is!!!!!
        cobra = trainSurrogates(cobra)  # side effect: constraintSurrogates, fitnessSurrogate

        ##########################################################
        # STEP6.2: Determine Distance requirement                #
        ##########################################################
        
        gama = cobra['XI'][((len(cobra['A'])-nRepair) % len(cobra['XI']))]
        
        ro = gama*cobra['l']
        
        ##########################################################
        # STEP6.3: Select Iterate                                #
        ##########################################################
        
        ptm = time.time()
        
        subMin = []
        
        verboseprint(cobra['verbose'], False, cobra['seqOptimizer']+" optimization on surrogate ...")
        
        if cobra['sac']['RS']:
            cobra = RandomStart(cobra)
            xStart = cobra['xStart']
        else:
            xStart = cobra['xStart']
            
        if cobra['seqOptimizer']=='COBYLA':            
            cons = []
            cons.append({'type':'ineq','fun':gCOBRA})
            
            for factor in range(len(cobra['lower'])):
                lower = cobra['lower'][factor]
                l = {'type':'ineq','fun': lambda x, lb=lower, i=factor: x[i]-lb}
                cons.append(l)
            for factor in range(len(cobra['upper'])):
                upper = cobra['upper'][factor]
                u = {'type':'ineq','fun': lambda x, ub=upper, i=factor: ub-x[i]}
                cons.append(u)
            
            opts = {'maxiter':cobra['seqFeval'], 'tol':cobra['seqTol']}
            subMin = optimize.minimize(subProb2, xStart, constraints=cons, options=opts, method='COBYLA')
            subMin['feval'] = subMin['nfev']
        else:
            raise ValueError("To Be implemented")
        
        cobra['optimizationTime'] = np.append(cobra['optimizationTime'], time.time()-ptm)
        verboseprint(cobra['verbose'], False, "finished "+str(subMin['feval'])+" iterations, "+str(cobra['optimizationTime'][len(cobra['optimizationTime'])-1])+" sec.")
        
        if penaF[0]*penaF[1] < penaF[2]:
            penaF[0] = penaF[0]*penaF[1]
            
        sigmaD = checkDistanceReq(subMin, fitnessSurrogate, sigmaD, CHECKDIST)
        cobra['progressCount'] = cobra['progressCount']+1
        cobra['fe'] = cobra['fe']+1
        
        xNew = subMin['x']
        xNew = np.maximum(xNew, cobra['lower'])
        xNew = np.minimum(xNew, cobra['upper'])
        newPredY = getPredY(xNew, fitnessSurrogate, cobra)
        
        if cobra['trueFuncForSurrogates']:
            newPredY = fn(xNew)[0]
        
        predY = np.append(predY, newPredY)
        predVal = np.append(predVal, subMin['fun'])
        feval = np.append(feval, subMin['feval'])
        optimizerConvergence = np.append(optimizerConvergence, subMin['status'])
        
        
        cobra['predC'] = np.vstack((cobra['predC'], getConstraintPrediction(xNew,constraintSurrogates)))
        
        xNewTemp = xNew
        xNewEval = fn(xNewTemp)
        
        if (((len(cobra['A'])%cobra['sac']['onlineFreqPLOG'])==0 and cobra['sac']['onlinePLOG']) 
            or (len(cobra['A'])==cobra['initDesPoints'])):
            
            newPredY1 = getPredY1(xNew, fitnessSurrogate1, cobra)
            newPredY2 = getPredY1(xNew, fitnessSurrogate2, cobra)
            newErr1 = np.abs(newPredY1-xNewEval[0])
            newErr2 = np.abs(plogReverse(newPredY2)-xNewEval[0])
            
            err1 = np.append(err1, newErr1)
            err2 = np.append(err2, newErr2)
            errRatio = err1/err2
            
            if np.isinf(newErr2):
                errRatio[len(errRatio)-1] = 0
            elif np.isinf(newErr1):
                errRatio[len(errRatio)-1] = np.inf
            
            cobra['pEffect'] = np.log10(np.nanpercentile(errRatio,75))
        
        newNumViol = np.sum(xNewEval[1:] > cobra['conTol'])
        feas = np.append(feas, newNumViol<1)
        newNumPred = np.sum(cobra['predC'][len(cobra['predC'])-1] > cobra['conTol'])
        feasPred = np.append(feasPred, newNumPred<1)
        
        if max(0,max(xNewEval[1:]))>cobra['conTol']:
            newMaxViol = max(0, max(xNewEval[1:]))
        else:
            newMaxViol = 0
        
        cobra = updateInfoAndCounters(cobra, xNew, xNewEval, newNumViol, newMaxViol, phase)
        
        Cfeas, EPS = adjustMargins(Cfeas, Tfeas, Cinfeas, Tinfeas, EPS, cobra['epsilonMax'])
        
        n = len(cobra['A'])
        
        cobra = updateSaveCobra(cobra, xNew, feas, feasPred, feval, optimizerConvergence, 
                                predY, predVal, subMin, sigmaD, penaF, gama, EPS)

    return(cobra)
