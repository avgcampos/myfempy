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
from colorama import Fore, Back, Style
    
#%% MYFEMPY DYNAMIC STEADY STATE SOLVE
def eig(usrlog,datamesh,F,freedof,KG,MG):
    if usrlog.solver_def == 'harnmodal':
        if usrlog.solver_opt == 'mode':
            modeStart = int(usrlog.solver_start)
            modeEnd = int(usrlog.solver_end)
            modeStep = int(usrlog.solver_step)
            
            if modeEnd == 0:
                print("input erro: solver_opt1 don't defined")
                modeEnd = 1

            if usrlog.mat_opt == 'lumped':
                W, V = sl.eigh((KG[:,freedof][freedof,:]).toarray(),(MG[:,freedof][freedof,:]).toarray())          
                modeEnd = len(W)
                U = sp.csc_matrix((datamesh[3],modeEnd))
                U[freedof,:] = V
            else:
                U = sp.csc_matrix((datamesh[3],modeEnd))
                W, U[freedof,:] = ssl.eigsh(KG[:,freedof][freedof,:],modeEnd,MG[:,freedof][freedof,:],sigma=1,which='LM')
            
            Wlist = np.arange(0,modeEnd+1)
            Wrad = np.sqrt(W)
            Whz = Wrad/(2*np.pi)
            U2 = np.concatenate((Wlist[1:, np.newaxis],Wrad[:, np.newaxis],Whz[:, np.newaxis]),axis=1)
            U3 = np.zeros((1,1))
        
        elif usrlog.solver_opt == 'fqrp':
            
            freqStart = (2*np.pi)*int(usrlog.solver_start)
            freqEnd = (2*np.pi)*int(usrlog.solver_end)
            freqStep = int(usrlog.solver_step)
            
            w_range =  np.linspace(freqStart,freqEnd,freqStep)
            U3 = sp.csc_matrix((datamesh[3],freqStep))
            
            # dofs = datamesh[3]
            # ith = np.zeros((freqStep))
            # jth = np.zeros((dofs))
            # val = np.zeros((len(freedof)))
            for ww in range(freqStep):
                Wn = w_range[ww]
                Dw = (KG[:,freedof][freedof,:])-(Wn**2)*(MG[:,freedof][freedof,:])
                U3[freedof,ww] = ssl.spsolve(Dw,F[freedof,:])
                # U_modal = ssl.spsolve(Dw,F[freedof,:])
            
            
            w_range = w_range/(2*np.pi)
            U2 = np.array([w_range])
            U = U3[:,-1]
     # elif usrlog.solver_def == 'supermode':
        # modeStart=int(usrlog.solver_start)
        # modeEnd=int(usrlog.solver_end)
        # modeStep=int(usrlog.solver_step)
                
    return U, U2, U3


    
#%% SOLVE SOLUTION
def solvesetup(usrlog,datamesh,coord,inci,tabmat,tabgeo,forcenodeaply,F,freedof):
    
    if usrlog.solver_typ == 'staticlinear':
        from myfempy.lib.myfempy_library import stifness_matrix
    
        start = time.time()
        KG = stifness_matrix(usrlog,datamesh,coord,inci,tabmat,tabgeo,forcenodeaply)
        end = time.time()
        kg_mem_size = sys.getsizeof(KG)/1e6
        assebtime = end - start
        print(Style.RESET_ALL)
        print(Fore.RED + Style.DIM+'\nKG: ',kg_mem_size,' Mb')
        print(Fore.RED + Style.DIM+'\nASSEMBLE TIME ',assebtime)
        
        start = time.time()
        U, U2, U3 = StaticLinear(usrlog,datamesh,F,freedof,KG)
        end = time.time()
        solvetime = end - start
        print(Style.RESET_ALL)
        print(Fore.RED + Style.DIM+'\nSOLVE TIME ',solvetime)
        
    # elif usrlog.solver_typ == 'staticnonlin':
    
    elif usrlog.solver_typ == 'dynsteadystt':
        
        from  myfempy.lib.myfempy_library import stifness_matrix, mass_matrix
    
        start = time.time()
        KG = stifness_matrix(usrlog,datamesh,coord,inci,tabmat,tabgeo,forcenodeaply)
        MG = mass_matrix(usrlog,datamesh,coord,inci,tabmat,tabgeo,forcenodeaply)
        end = time.time()
        assebtime = end - start
        kg_mem_size = sys.getsizeof(KG)/1e6
        print(Style.RESET_ALL)
        print(Fore.RED + Style.DIM+'\nSIZE KG: ',kg_mem_size,' Mb')
        print(Fore.RED + Style.DIM+'\nASSEMBLY TIME ',assebtime)
         
        start = time.time()
        U, U2, U3 = DynSteadyStt(usrlog,datamesh,F,freedof,KG,MG)
        end = time.time()
        solvetime = end - start
        print(Style.RESET_ALL)
        print(Fore.RED + Style.DIM+'\nSOLVE TIME ', solvetime)
            
    time_list = np.array([assebtime,solvetime,kg_mem_size])
    return KG, U, U2, U3, time_list