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
import scipy.sparse.linalg as spla
import scipy.linalg as sl
import time

# %% MYFEMPY STATIC SOLVE


def EIG(fulldofs, stiffness, mass, forcelist, freedof, solverset):

    plotset = dict()
    postprocset = dict()
    U = np.zeros((fulldofs, solverset['end']))

    modeEnd = solverset['end']

    startstep = time.time()
    W, U[freedof, :] = spla.eigsh(stiffness[:, freedof][freedof, :],
                                 modeEnd, mass[:, freedof][freedof, :], sigma=1, which='LM')
    endstep = time.time()
    print('\nSTEP --: SUCCESSFUL CONVERGED\n')

    print('\ TIME SPEND: ', endstep-startstep, ' SEC')

    Wlist = np.arange(0, modeEnd+1)
    Wrad = np.sqrt(W)
    Whz = Wrad/(2*np.pi)
    w_range = np.concatenate(
        (Wlist[1:, np.newaxis], Wrad[:, np.newaxis], Whz[:, np.newaxis]), axis=1)

    return U, w_range


def FRF(fulldofs, stiffness, mass, forcelist, freedof, solverset):

    twopi = 2*np.pi
    freqStart = (twopi)*solverset['start']
    freqEnd = (twopi)*solverset['end']
    freqStep = solverset['nsteps']

    U = np.zeros((fulldofs, solverset['nsteps']))
    w_range = np.linspace(freqStart, freqEnd, freqStep)

    startstep = time.time()
    for ww in range(freqStep):
        Wn = w_range[ww]
        
        Dw = (stiffness[:, freedof][freedof, :]) - \
            (Wn**2)*(mass[:, freedof][freedof, :])

        U[freedof, ww] = spla.spsolve(A=Dw, b=forcelist[freedof,:])
        
        # U[freedof, ww], info = spla.bicgstab(A=Dw, b=forcelist[freedof, :].toarray(), tol=solverset['TOL'])
        
        # if info > 0:
        #     pass
        # elif info < 0:
        #     print('ILLEGAL INPUT OR BREAKDOWN')
        # else:
        #     pass
        
    endstep = time.time()



    print('\ TIME SPEND: ', endstep-startstep, ' SEC')

    w_range = w_range/(twopi)

    return U, w_range
