# -*- coding: utf-8 -*-
"""
Created on Fri Nov 24 10:37:59 2017

@author: r.dewinter
"""

def defaultSAC(DOSAC=1):
    '''
    Default settings for the SACOBRA part of COBRA.

    Sets suitable defaults for the SACOBRA part of COBRA. 
    With the call setOpts(mySAC,defaultSAC()) it is possible to extend a partial list 
    mySAC to a list containing all sac-elements (the missing ones are taken from 
    defaultSAC()).
    
    For backward compatibility, a logical DOSAC (deprecated) is mapped from FALSE to 0 and 
    from TRUE to 1.
    
     @param DOSAC  [0|1|2] 0: COBRA-R settings (turn off SACOBRA), 1: SACOBRA settings, 2: SACOBRA settings 
                   with fewer parameters and more online adujustement (aFF and aCF are done parameter free).
                           
     @return a list with the following elements (the values in parantheses [|] are the values 
       for DOSAC=[0|1|2]): 
         RS   flag for random start algorithm [FALSE|TRUE|TRUE]  
         RStype type of the function to calculate probability to start the internal optimizer 
                     with a random starting point[NA|SIGMOID|CONSTANT] (see function 
                     RandomStart in SACOBRA.R)
         RS max maximum probability of a random start when RStype is selected as SIGMOID (see RandomStart in SACOBRA.R) 
         if RStype is CONSTANT then random start is done with a constant probability determined from mean(c(RSmax,RSmin)) [NA|0.3|0.3]
         RSmin manimum probability of a random start when RStype is selected as SIGMOID (see RandomStart in SACOBRA.R)[NA|0.05|0.05] 
         aDRC   flag for automatic DRC adjustment [FALSE|TRUE|TRUE]  
         aFF   flag for automatic objective function transformation [FALSE|TRUE|TRUE]  
         aCF   flag for automatic constraint function transformation [FALSE|TRUE|TRUE]  
         TFRange   threshold, if FRange is larger, then apply automatic objective
             function transformation (see plog). [Inf|1e+05|-1]  
         TGR   threshold, if GRatio is larger, then apply automatic constraint
             function transformation. GRatio is the ratio "largest GRange /
             smallest GRange" where GRange is the min-max range of a specific constraint.
             If TGR < 1, then the transformation is always performed.  [Inf|1e+03|-1].
         CS   threshold unsuccessful iterations for random start algorithm. If  CS
             iterations in a row do not improve the ever-best feasible solution, then perform a 
             restart. [10|10|10]  
         adaptivePLOG   (experimental) flag for objective function transformation with plog, 
             where the parameter pShift is adapted during iterations. 
             [FALSE|FALSE|FALSE] 
         onlinePLOG  flag for online decision making wether use plog or not according to p-effect plog. 
             [FALSE|FALSE|TRUE] 
         pEffectInit Initial pEffect value when using onlinePLOG. If pEffectInit >= 2 then the initial model is built after plog transformation.
             [NA|NA|2]
       
     @seealso   cobraInit, cobraPhaseII
     @author Samineh Bagheri, Cologne Univeristy of Applied Sciences
    '''
    if DOSAC==1:
        sac = dict()
        sac['RS'] = True
        sac['RStype'] = "SIGMOID"  
        sac['RSAUTO'] = False
        sac['RSmax'] = 0.3  #maximum probability of a random start
        sac['RSmin'] = 0.05   #minimum probability of a random start
        sac['aDRC'] = True
        sac['aFF'] = True
        sac['TFRange'] = 1e+05
        sac['aCF'] = True
        sac['TGR'] = 1e+03     #GR threshold to perform constraint modification,
#                                 if TGR < 1 transformation is always performed,
#                                 Inf value means transfomration is never performed
#                                 a positive value laarger than 1 means: 
#                                 the transfomartion is only performed for problems with a GR larger than TGR
        sac['conPLOG'] = False  #perform plog transformation for all constraints (In testing phase) 
        sac['conFitPLOG'] = False #perform plog transformation for all constraints and fitness (In testing phase)
        sac['adaptivePLOG'] = False #the effectivity of adaptivePLOG is not fully proved therefore I keep it as FALSE for now
        sac['onlinePLOG'] = False
        sac['onlineFreqPLOG'] = 10 # number of iterations in a row which after that the online DOPLOG check is done
        sac['pEffectInit'] = 0
        sac['minMaxNormal'] = False
        sac['onlineMinMax'] = False
        sac['Cs'] = 10
    else:
        raise ValueError('not implemented yet')
    return(sac)
    
def setOpts(opts, defaultOpt):
    for key in opts:
        if key in defaultOpt:
            opts[key] = defaultOpt[key]
    return opts