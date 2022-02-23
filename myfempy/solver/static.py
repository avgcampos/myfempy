# -*- coding: utf-8 -*-
"""
Created on Thu Mar 25 15:33:44 2021
@author: ANTONIO VINICIUS GARCIA CAMPOS
@version: v20
_______________________________________________________________________________
 ~~~~~~~~~~                  SOLVER LINEAR ELASTICO                  ~~~~~~~~~~

> FUNCIONALIDADES
--- ENTRADAS:
--- SAIDA:  

===============================================================================

> ATUALIZACOES DA VERSAO:
--- 
--- 
_______________________________________________________________________________
"""
import sys
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as ssl
import scipy.linalg as sl
import time
from colorama import Fore, Back, Style

from .setup import step_setting
from ..io.miscel import loading_bar_v1

#%% MYFEMPY STATIC SOLVE
def solve(datamesh,forcelist,freedof,KG,solverset):  
    U = np.zeros((datamesh['fulldof'],1))
    if solverset['type'] == 'direct':
        U[freedof,0] = ssl.spsolve(KG[:,freedof][freedof,:],forcelist[freedof,solverset['step']])
    
    elif solverset['type'] == 'iterative':
        U[freedof,0], info = ssl.bicg(KG[:,freedof][freedof,:],forcelist[freedof,solverset['step']],tol=solverset['tol'])
        
        print('\n')
        if info > 0:
            print('Step --',solverset['step'],': Not converged to tolerance achieved, exceeded number of iterations...\n')
        else:
            print('Step --',solverset['step'],': Successful converged...\n')
    print('>> Solver exit -- Solution finished...\n')
    return U

   
#%% SOLVE SOLUTION
def alglinear(usrlog,datamesh,coord,inci,tabmat,tabgeo,forcenodeaply,forcelist,freedof):
    from myfempy.lib.myfempy_library import stifness_matrix

    start = time.time()
    KG = stifness_matrix(usrlog,datamesh,coord,inci,tabmat,tabgeo,forcenodeaply)
    end = time.time()
    kg_mem_size = sys.getsizeof(KG)/1e6
    assebtime = end - start
    print(Style.RESET_ALL)
    print(Fore.RED + Style.DIM+'\nKG: ',kg_mem_size,' Mb')
    print(Fore.RED + Style.DIM+'\nASSEMBLE TIME ',assebtime)
    
        
    stepstart = int(usrlog.solver_start)
    stepend = int(usrlog.solver_end)
    stepstep = int(usrlog.solver_step)
    
    stepset = step_setting(usrlog, stepstart, stepend, stepstep)
    U0 = np.zeros((datamesh[3],1))
    U1 = np.zeros((datamesh[3],1))
    U = sp.csc_matrix((datamesh[3],stepset))
    
    start = time.time()
    for ustep in range(stepset):
        loading_bar_v1(100*((ustep+1)/stepset),'SOSTEP')
        solvecfg = {'type' : usrlog.solver_def,
                    'step' : ustep}
        
        U1 = solve(datamesh, forcelist, freedof, KG, solvecfg)
        U1 += U0
        U[:, ustep] = U1
        U0 = U1
        
    U2 = np.zeros((stepset,1))
    U3 = np.zeros((stepset,1))      
    end = time.time()
    solvetime = end - start
    print(Style.RESET_ALL)
    print(Fore.RED + Style.DIM+'\nSOLVE TIME ',solvetime)
                    
    time_list = np.array([assebtime,solvetime,kg_mem_size])
    return KG, U, U2, U3, time_list

#%%
#def nonlinear():