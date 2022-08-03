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

from myfempy.core.assembler import assembler, loads
from myfempy.core.solverset import get_constrains_dofs, step_setting
from myfempy.tools.tools import print_console


def get_solve(solver_type):

    if solver_type == 'SLD':
        from myfempy.core.staticlinear import SLD
        return SLD

    elif solver_type == 'SLI':
        from myfempy.core.staticlinear import SLI
        return SLI

    elif solver_type == 'SLIPRE':
        from myfempy.core.staticlinear import SLIPRE
        return SLIPRE

    elif solver_type == 'EIG':
        from myfempy.core.vibralinear import EIG
        return EIG

    elif solver_type == 'FRF':
        from myfempy.core.vibralinear import FRF
        return FRF

    else:
        print('Erro Import Solver')


# %% SOLVE SOLUTION
def gen_static_solution(solverset, modelinfo):

    print_console('solver')

    solve = get_solve(solverset['SOLVER'])

    start = time.time()
    KG = assembler(modelinfo, key='stiffness')
    end = time.time()
    kg_mem_size = 0  # sys.getsizeof(KG.toarray())/1e6
    assembly_time = end - start
    print(' ')
    print('STIFFNESS SIZE: ', kg_mem_size, ' MB')
    print('ASSEMBLY FULL TIME SPEND ', assembly_time, ' SEC')

    F, KG = loads(modelinfo, KG)
    freedof, fixedof = get_constrains_dofs(modelinfo)

    solverset['nsteps'] = step_setting(solverset['STEPSET'])
    solverset['coord'] = modelinfo['coord']
    # solverset["TRACKER"] =

    fulldofs = (modelinfo["nodedof"][0])*len(modelinfo["coord"])
    start = time.time()
    U = solve(fulldofs, KG, F, freedof, solverset)
    end = time.time()
    time_spend = end - start
    print('\nSOLVE FULL TIME SPEND ', time_spend, ' SEC')

    solvestatus = {'timeasb': assembly_time,
                   'timesim': time_spend,
                   'kgsize': kg_mem_size}

    solution = {'U': U,
                'solvestatus': solvestatus}

    print_console('thank')

    return solution


def gen_dynamic_solution(solverset, modelinfo):

    print_console('solver')

    solve = get_solve(solverset['SOLVER'])

    start = time.time()
    KG = assembler(modelinfo, key='stiffness')
    MG = assembler(modelinfo, key='mass')
    end = time.time()
    kg_mem_size = 0  # sys.getsizeof(KG.toarray())/1e6
    assembly_time = end - start
    print(' ')
    print('STIFFNESS SIZE: ', kg_mem_size, ' MB')
    print('ASSEMBLY FULL TIME SPEND ', assembly_time, ' SEC')

    F, KG = loads(modelinfo, KG)
    freedof, fixedof = get_constrains_dofs(modelinfo)

    solverset['start'] = solverset['STEPSET']['start']
    solverset['end'] = solverset['STEPSET']['end']
    solverset['nsteps'] = step_setting(solverset['STEPSET'])
    solverset['coord'] = modelinfo['coord']
    # solverset["TRACKER"] =

    fulldofs = (modelinfo["nodedof"][0])*len(modelinfo["coord"])
    start = time.time()
    U, w_range = solve(fulldofs, KG, MG, F, freedof, solverset)
    end = time.time()
    time_spend = end - start
    print('\nSOLVE FULL TIME SPEND ', time_spend, ' SEC')

    solvestatus = {'timeasb': assembly_time,
                   'timesim': time_spend,
                   'kgsize': kg_mem_size}

    solution = {'U': U,
                'FREQ': w_range,
                'solvestatus': solvestatus,
                }

    print_console('thank')

    return solution
