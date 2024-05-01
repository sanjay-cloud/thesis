# -*- coding: utf-8 -*-
"""
Created on Fri Nov 24 09:56:40 2017

@author: r.dewinter
"""

def defaultRI(repairMargin=1e-2):
    '''
    defaultRI
    
    Default settings for repairInfeasRI2 and repairChootinan.
    
    Sets suitable defaults for the repair-infeasible part of COBRA. cr
    With the call setOpts(myRI,defaultRI()) it is possible to extend a partial list 
    myRI to a list containing all ri-elements (the missing ones are taken from 
    defaultRI())  
    
    Detail:
    A solution x is said to be epsilon-feasible for constraint function f, if 
          f(x)+epsilon leq 0 
          
    The infeasibility of a solution is its maximum constraint violation 
    (0 for a feasible solution).     
    
     @param repairMargin  [1e-2] repair only solutions whose infeasibility is less than this margin 
    
     @return a list with the following elements:
       
         RIMODE{  [2] one out of 0,1,2,3  with 0,1: deprecated older versions of RI2, 
             2: the recommended RI2-case, see repairInfeasRI2, 
             3: Chootinan's method, see repairChootinan  
         eps1  [1e-4] include all constraints not eps1-feasible into the repair mechanism  
         eps2  [1e-4] selects the solution with the shortest shift among all random 
             realizations which are eps2-feasible  
         q  [3.0] draw coefficients alpha_k from uniform distribution U[0,q]  
         mmax  [1000] draw mmax random realizations  
         repairMargin  repair only solutions whose infeasibility is less than this margin. 
         repairOnlyFresBetter{  [FALSE] if TRUE, then repair only iterates with cr
             fitness < so-far-best-fitness + marFres  
         marFres  [0.0] only relevant if repairOnlyFresBetter==TRUE 
       
         
     @seealso   repairInfeasRI2, repairChootinan
    
     @author Wolfgang Konen, Cologne Univeristy of Applied Sciences
    '''
    ri = dict()
    ri['RIMODE'] = 2 # 0: OLD RI, 1: RI w/o epsilon-feasibility, 2: RI2, the recommended case 3: repairChootinan
    ri['eps1'] = 1e-4 # include all constraints not eps1-feasible into the repair mechanism
    ri['eps2'] = 1e-4 # selectBest() selects the solution with the shortest shift among all random realizations which are eps2-feasible 
    ri['q'] = 3.0 # draw alpha_k from uniform distribution U[0,q]
    ri['mmax'] = 1000 # draw mmax random realizations
    ri['OLD'] = False # TRUE: activate the old repairInfeasible (before 2014-09-29)
    ri['kappa'] = 1.2 # (only OLD) if =1.0: try to step directly to the true boundary, if >1.0: move q bit further into the feasible region
    ri['repairMargin'] = repairMargin # repair only solutions whose infeasibility is less than this margin
    ri['repairOnlyFresBetter'] = False # if repairOnlyFresBetter=TRUE, then repair only iterates with fitness < so-far-best-fitness + marFres
    ri['marFres'] = 0 # only relevant if repairOnlyFresBetter==TRUE 
    return(ri)