# -*- coding: utf-8 -*-
"""
Created on Fri Nov 24 10:29:51 2017

@author: r.dewinter
"""

def defaultTR():
    '''
    defaultTR.r
    Default settings for the trust-region part of COBRA.
     
    Sets default values for the trust-region part cobra$TRlist of COBRA.  
    With the call setOpts(myTR,defaultTR()) it is possible to extend a partial list 
    myTR to a list containing all TR-elements (the missing ones are taken from 
    defaultTR()).
    
    
    @return a list with the following elements 
        radiMin Minimum fraction of the width of the search space to be used as radius of the trust region [c(0:1)]
        radiMax Maximum fraction of the width of the search space to be used as radius of the trust region [c(0:1)]
        radiInit Initial radius of trust region
    '''
    tr = dict()
    tr['shape'] = "cube"
    tr['radiMin'] = 0.01
    tr['radiMax'] = 0.8
    tr['radiInit'] = 0.1
    tr['center'] = "xbest"
    return(tr)