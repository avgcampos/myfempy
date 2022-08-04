#!/usr/bin/env python
from myfempy.plots.plotxy import tracker_plot
from myfempy.tools.tools import loading_bar_v1
import time
import scipy.sparse.linalg as spla
import numpy as np
__doc__ = """
Static Linear Solver
"""


def sld(fulldofs, stiffness, forcelist, freedof, solverset):
    plotset = dict()
    postprocset = dict()
    U0 = np.zeros((fulldofs, 1))
    U1 = np.zeros((fulldofs, 1))
    U = np.zeros((fulldofs, solverset['nsteps']))
    loading_bar_v1(0, 'SOLVER')
    for step in range(solverset['nsteps']):
        loading_bar_v1(100*((step+1)/solverset['nsteps']), 'SOLVER')
        startstep = time.time()
        U1[freedof, 0] = spla.spsolve(
            A=stiffness[:, freedof][freedof, :], b=forcelist[freedof, step])
        endstep = time.time()
        print('\nSTEP --', step, ': SUCCESSFUL CONVERGED\n')
        print('\nSOLVE STEP '+str(step),
              '\nTIME SPEND: ', endstep-startstep, ' SEC\n')
        U1[freedof, 0] += U0[freedof, 0]
        U[freedof, step] = U1[freedof, 0]
        U0[freedof, 0] = U1[freedof, 0]
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


def sli(fulldofs, stiffness, forcelist, freedof, solverset):
    plotset = dict()
    postprocset = dict()
    U0 = np.zeros((fulldofs, 1))
    U1 = np.zeros((fulldofs, 1))
    U = np.zeros((fulldofs, solverset['nsteps']))
    loading_bar_v1(0, 'SOLVER')
    for step in range(solverset['nsteps']):
        loading_bar_v1(100*((step+1)/solverset['nsteps']), 'SOLVER')
        startstep = time.time()
        U1[freedof, 0], info = spla.bicgstab(A=stiffness[:, freedof][freedof, :], b=forcelist[freedof, step].toarray(),
                                             x0=U0[freedof, 0], tol=solverset['TOL'])
        endstep = time.time()
        if info > 0:
            print('\nSTEP --', step, ': CONVERGED TO TOLERANCE NOT ACHIEVED\n')
        elif info < 0:
            print('\nSTEP --', step, ': ILLEGAL INPUT OR BREAKDOWN\n')
        else:
            print('\nSTEP --', step, ': SUCCESSFUL CONVERGED\n')

        print('\nSOLVE STEP '+str(step),
              '\nTIME SPEND: ', endstep-startstep, ' SEC\n')
        U1[freedof, 0] += U0[freedof, 0]
        U[freedof, step] = U1[freedof, 0]
        U0[freedof, 0] = U1[freedof, 0]
        if "TRACKER" in solverset.keys():
            if solverset["TRACKER"]['show'] == True:
                plotset['step'] = step+1
                plotset['val_list'] = U1
                plotset['fignumb'] = 1
                postprocset["TRACKER"] = solverset["TRACKER"]
                coord = solverset['coord']
                tracker_plot(postprocset, plotset, coord)
    return U


def slipre(fulldofs, stiffness, forcelist, freedof, solverset):
    plotset = dict()
    postprocset = dict()
    U0 = np.zeros((fulldofs, 1))
    U1 = np.zeros((fulldofs, 1))
    U = np.zeros((fulldofs, solverset['nsteps']))
    print('PRECONDITIONING M MATRIX\n')
    Nshape = np.size(freedof)
    sA = stiffness[:, freedof][freedof, :]
    def A_ilu(x): return spla.spsolve(sA, x)
    M = spla.LinearOperator((Nshape, Nshape), A_ilu)
    loading_bar_v1(0, 'SOLVER')
    for step in range(solverset['nsteps']):
        loading_bar_v1(100*((step+1)/solverset['nsteps']), 'SOLVER')
        startstep = time.time()
        U1[freedof, 0], info = spla.gmres(A=sA, b=forcelist[freedof, step].toarray(),
                                          x0=U0[freedof, 0], tol=solverset['TOL'], M=M)
        endstep = time.time()
        if info > 0:
            print('\nSTEP --', step, ': CONVERGED TO TOLERANCE NOT ACHIEVED\n')
        elif info < 0:
            print('\nSTEP --', step, ': ILLEGAL INPUT OR BREAKDOWN\n')
        else:
            print('\nSTEP --', step, ': SUCCESSFUL CONVERGED\n')
        print('\nSOLVE STEP '+str(step),
              '\nTIME SPEND: ', endstep-startstep, ' SEC\n')
        U1[freedof, 0] += U0[freedof, 0]
        U[freedof, step] = U1[freedof, 0]
        U0[freedof, 0] = U1[freedof, 0]
        if "TRACKER" in solverset.keys():
            if solverset["TRACKER"]['show'] == True:
                plotset['step'] = step+1
                plotset['val_list'] = U1
                plotset['fignumb'] = 1
                postprocset["TRACKER"] = solverset["TRACKER"]
                coord = solverset['coord']
                tracker_plot(postprocset, plotset, coord)
    return U
