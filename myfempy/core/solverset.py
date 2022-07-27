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

import numpy as np

#------------------------------------------------------------------------------
#%% vetor cond. contorno
def get_constrains_dofs(modelinfo):
    # numero de tipos de restriçoes fixas
        
    ntbc = modelinfo['constrains'].shape[0] 

    # dofs fixos do sistema
    fixedof = np.zeros((1,modelinfo["nodedof"][0]*len(modelinfo["coord"])))
    
    if modelinfo['nodedof'][0] == 2:
        if (modelinfo['elemid'][0] == 110) or (modelinfo['elemid'][0] == 120) or (modelinfo['elemid'][0] == 210) or (modelinfo['elemid'][0] == 220):
            for ii in range(ntbc):
                no = int(modelinfo['constrains'][ii,1])
                if int(modelinfo['constrains'][ii,0]) == 0:
                    fixedof[0,modelinfo['nodedof'][0]*no-2] = modelinfo['nodedof'][0]*no-1
                    fixedof[0,modelinfo['nodedof'][0]*no-1] = modelinfo['nodedof'][0]*no
                elif int(modelinfo['constrains'][ii,0]) == 1:
                    fixedof[0,modelinfo['nodedof'][0]*no-2] = modelinfo['nodedof'][0]*no-1
                elif int(modelinfo['constrains'][ii,0]) == 2:
                    fixedof[0,modelinfo['nodedof'][0]*no-1] = modelinfo['nodedof'][0]*no

        elif modelinfo['elemid'][0] == 130:
            for ii in range(ntbc):
                no = int(modelinfo['constrains'][ii,1])
                if int(modelinfo['constrains'][ii,0]) == 0:
                    fixedof[0,modelinfo['nodedof'][0]*no-2] = modelinfo['nodedof'][0]*no-1
                    fixedof[0,modelinfo['nodedof'][0]*no-1] = modelinfo['nodedof'][0]*no
                elif int(modelinfo['constrains'][ii,0]) == 2:
                    fixedof[0,modelinfo['nodedof'][0]*no-2] = modelinfo['nodedof'][0]*no-1
                elif int(modelinfo['constrains'][ii,0]) == 6:
                    fixedof[0,modelinfo['nodedof'][0]*no-1] = modelinfo['nodedof'][0]*no             
    

    elif modelinfo['nodedof'][0] == 3:
        if modelinfo['elemid'][0] == 140:
            for ii in range(ntbc):
                no = int(modelinfo['constrains'][ii,1])
                if int(modelinfo['constrains'][ii,0]) == 0:
                    fixedof[0,modelinfo['nodedof'][0]*no-3] = modelinfo['nodedof'][0]*no-2
                    fixedof[0,modelinfo['nodedof'][0]*no-2] = modelinfo['nodedof'][0]*no-1
                    fixedof[0,modelinfo['nodedof'][0]*no-1] = modelinfo['nodedof'][0]*no
                elif int(modelinfo['constrains'][ii,0]) == 1:
                    fixedof[0,modelinfo['nodedof'][0]*no-3] = modelinfo['nodedof'][0]*no-2
                elif int(modelinfo['constrains'][ii,0]) == 2:
                    fixedof[0,modelinfo['nodedof'][0]*no-2] = modelinfo['nodedof'][0]*no-1
                elif int(modelinfo['constrains'][ii,0]) == 6:
                    fixedof[0,modelinfo['nodedof'][0]*no-1] = modelinfo['nodedof'][0]*no  
        
        elif (modelinfo['elemid'][0] == 310) or (modelinfo['elemid'][0] == 320):
            for ii in range(ntbc):
                no = int(modelinfo['constrains'][ii,1])
                if int(modelinfo['constrains'][ii,0]) == 0:
                    fixedof[0,modelinfo['nodedof'][0]*no-3] = modelinfo['nodedof'][0]*no-2
                    fixedof[0,modelinfo['nodedof'][0]*no-2] = modelinfo['nodedof'][0]*no-1
                    fixedof[0,modelinfo['nodedof'][0]*no-1] = modelinfo['nodedof'][0]*no
                elif int(modelinfo['constrains'][ii,0]) == 1:
                    fixedof[0,modelinfo['nodedof'][0]*no-3] = modelinfo['nodedof'][0]*no-2
                elif int(modelinfo['constrains'][ii,0]) == 2:
                    fixedof[0,modelinfo['nodedof'][0]*no-2] = modelinfo['nodedof'][0]*no-1
                elif int(modelinfo['constrains'][ii,0]) == 3:
                    fixedof[0,modelinfo['nodedof'][0]*no-1] = modelinfo['nodedof'][0]*no
    
    elif modelinfo['nodedof'][0] == 6:
        if modelinfo['elemid'][0] == 141:
            for ii in range(ntbc):
                no = int(modelinfo['constrains'][ii,1])
                if int(modelinfo['constrains'][ii,0]) == 0:
                    fixedof[0,modelinfo['nodedof'][0]*no-6] = modelinfo['nodedof'][0]*no-5
                    fixedof[0,modelinfo['nodedof'][0]*no-5] = modelinfo['nodedof'][0]*no-4
                    fixedof[0,modelinfo['nodedof'][0]*no-4] = modelinfo['nodedof'][0]*no-3
                    fixedof[0,modelinfo['nodedof'][0]*no-3] = modelinfo['nodedof'][0]*no-2
                    fixedof[0,modelinfo['nodedof'][0]*no-2] = modelinfo['nodedof'][0]*no-1
                    fixedof[0,modelinfo['nodedof'][0]*no-1] = modelinfo['nodedof'][0]*no
                elif int(modelinfo['constrains'][ii,0]) == 1:
                    fixedof[0,modelinfo['nodedof'][0]*no-6] = modelinfo['nodedof'][0]*no-5
                elif int(modelinfo['constrains'][ii,0]) == 2:
                    fixedof[0,modelinfo['nodedof'][0]*no-5] = modelinfo['nodedof'][0]*no-4
                elif int(modelinfo['constrains'][ii,0]) == 3:
                    fixedof[0,modelinfo['nodedof'][0]*no-4] = modelinfo['nodedof'][0]*no-3
                elif int(modelinfo['constrains'][ii,0]) == 4:
                    fixedof[0,modelinfo['nodedof'][0]*no-3] = modelinfo['nodedof'][0]*no-2
                elif int(modelinfo['constrains'][ii,0]) == 5:
                    fixedof[0,modelinfo['nodedof'][0]*no-2] = modelinfo['nodedof'][0]*no-1
                elif int(modelinfo['constrains'][ii,0]) == 6:
                    fixedof[0,modelinfo['nodedof'][0]*no-1] = modelinfo['nodedof'][0]*no
            
    
    fixedof = fixedof[np.nonzero(fixedof)]
    fixedof = fixedof - np.ones_like(fixedof)
    # fixedofs = np.unique(fixedofs)
    # todos os dofs do sistema
    alldof = np.arange(0,modelinfo["nodedof"][0]*len(modelinfo["coord"]),1,int)
    # dofs livres do sistema
    freedof = np.setdiff1d(alldof,fixedof)
    return freedof, fixedof



#------------------------------------------------------------------------------
#%% step setting
def step_setting(steps):
    
    start = steps['start']
    end = steps['end']
    substep = steps['step']
    
    if (end-start)==0:
        nsteps = int(end)
    else:
        nsteps = int((end-start)/substep)
    return nsteps


