# -*- coding: utf-8 -*-
"""
========================================================================
~~~ MODULO DE SIMULACAO ESTRUTURAL PELO METODO DOS ELEMENTOS FINITOS ~~~
       	                    __                                
       	 _ __ ___   _   _  / _|  ___  _ __ ___   _ __   _   _ 
       	| '_ ` _ \ | | | || |_  / _ \| '_ ` _ \ | '_ \ | | | |
       	| | | | | || |_| ||  _||  __/| | | | | || |_) || |_| |
       	|_| |_| |_| \__, ||_|   \___||_| |_| |_|| .__/  \__, |
       	            |___/                       |_|     |___/ 

~~~      Mechanical studY with Finite Element Method in PYthon       ~~~
~~~                PROGRAMA DE ANÁLISE COMPUTACIONAL                 ~~~
~~~              copyright @ 2022, all rights reserved               ~~~
========================================================================
"""
import sys
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as ssl
import scipy.linalg as sl
import time

from myfempy.tools.tools import loading_bar_v1
from myfempy.plots.plotxy import tracker_plot

#%% MYFEMPY STATIC SOLVE
def SLD(fulldofs,stiffness,forcelist,freedof,solverset):
    
    plotset = dict()
    postprocset = dict()
    U0 = np.zeros((fulldofs,1))
    U1 = np.zeros((fulldofs,1))
    U = np.zeros((fulldofs,solverset['nsteps']))
     
    loading_bar_v1(0,'SOLVER')
    for step in range(solverset['nsteps']):
        loading_bar_v1(100*((step+1)/solverset['nsteps']),'SOLVER')
        startstep = time.time()
        U1[freedof,0] = ssl.spsolve(stiffness[:,freedof][freedof,:], forcelist[freedof,step])
        endstep = time.time()
        print('\nSTEP --',step,': SUCCESSFUL CONVERGED\n')
       
        print('\nSOLVE STEP '+str(step),'\ TIME SPEND: ',endstep-startstep,' SEC')
        U1[freedof,0] += U0[freedof,0]
        U[freedof,step] = U1[freedof,0]
        U0[freedof,0] = U1[freedof,0]
        
        if "TRACKER" in solverset.keys():                
            if solverset["TRACKER"]['show'] == True:
               # for st in range(len(plotset[postprocset["TRACKER"]['result2plot']])):
                plotset['step'] = step+1
                plotset['val_list'] = U1
                plotset['fignumb'] = 1
                postprocset["TRACKER"] = solverset["TRACKER"]
                coord = solverset['coord']
                tracker_plot(postprocset, plotset, coord)
    
    return U



def SLI(fulldofs,stiffness,forcelist,freedof,solverset):  
    
    plotset = dict()
    postprocset = dict()
    U0 = np.zeros((fulldofs,1))
    U1 = np.zeros((fulldofs,1))
    U = np.zeros((fulldofs,solverset['nsteps']))
     
    loading_bar_v1(0,'SOLVER') 
    for step in range(solverset['nsteps']):
        loading_bar_v1(100*((step+1)/solverset['nsteps']),'SOLVER')
        startstep = time.time()
        U1[freedof,0], info = ssl.bicg(stiffness[:,freedof][freedof,:], forcelist[freedof,step], U0[freedof,0], solverset['TOL'])
        endstep = time.time()
        
        
        if info > 0:
            print('\nSTEP --',step,': CONVERGED TO TOLERANCE NOT ACHIEVED\n')
        elif info < 0:
            print('\nSTEP --',step,': ILLEGAL INPUT OR BREAKDOWN\n') 
        else:
            print('\nSTEP --',step,': SUCCESSFUL CONVERGED\n')
            
        print('\nSOLVE STEP '+str(step),'\ TIME SPEND: ',endstep-startstep,' SEC')    
        U1[freedof,0] += U0[freedof,0]
        U[freedof,step] = U1[freedof,0]
        U0[freedof,0] = U1[freedof,0]
        
        if "TRACKER" in solverset.keys():                
            if solverset["TRACKER"]['show'] == True:
                plotset['step'] = step+1
                plotset['val_list'] = U1
                plotset['fignumb'] = 1
                postprocset["TRACKER"] = solverset["TRACKER"]
                coord = solverset['coord']
                tracker_plot(postprocset, plotset, coord)
               
    return U
