from __future__ import annotations

import time

import numpy as np
# import jax.numpy as np
import scipy.sparse.linalg as spla

from myfempy.core.alglin import linsolve_gmres, linsolve_direct
from myfempy.core.solver import Solver


class StaticLinear(Solver):
    '''Static Linear Solver Class <ConcreteClassService>'''
    
    def getMatrixAssembler(Model, inci, coord, tabmat, tabgeo, intgauss):
                
        matrix = dict()
        startstep = time.time()
        matrix['stiffness'] = Solver.getMatrixAssembler(Model, inci, coord, tabmat, tabgeo, intgauss,  type_assembler = 'linear_stiffness')
        endstep = time.time()
        print("\nGLOBAL ASSEMBLY TIME ", "\ TIME SPEND: ", endstep - startstep, " SEC")
        return matrix

    def getLoadAssembler(loadaply, nodetot, nodedof):
        return Solver.getLoadAssembler(loadaply, nodetot, nodedof)
    
    def addMatrix(A, A_add, loc):
        return Solver.addMatrix(A, A_add, loc)
    
    def getConstrains(constrains, nodetot, nodedof):
        return Solver.getConstrains(constrains, nodetot, nodedof)
    
    def setSteps(steps):
        return Solver.setSteps(steps)
        

    def Solve(fulldofs, assembly, forcelist, freedof, solverset):
        """scipy generalized minimal residual iteration solver"""
        
        stiffness = assembly['stiffness']
        
        nsteps = StaticLinear.setSteps(solverset["STEPSET"])
        
        # plotset = dict()
        # postprocset = dict()
        U0 = np.zeros((fulldofs, 1))
        U1 = np.zeros((fulldofs, 1))
        U = np.zeros((fulldofs, nsteps))
        # print("PRECONDITIONING M MATRIX\n")
        Nshape = np.size(freedof)
        sA = stiffness[:, freedof][freedof, :]

        def A_ilu(x):
            return spla.spsolve(sA, x)

        M = spla.LinearOperator((Nshape, Nshape), A_ilu)
        # loading_bar_v1(0, "SOLVER")

        solution = dict()

        for step in range(nsteps):
            # loading_bar_v1(100 * ((step + 1) / solverset["nsteps"]), "SOLVER")
            startstep = time.time()
            
            # U1[freedof, 0], info = linsolve_gmres(sA, forcelist[freedof, step].toarray(), U0[freedof, 0], M)

            U1[freedof, 0] = linsolve_direct(stiffness[:, freedof][freedof, :], forcelist[freedof, step])
            
            # x, info = linsolve_Jaxgmres(sA, forcelist[freedof, step].toarray(), U0[freedof, 0], solverset["TOL"], M)            
            # U1.at[freedof, 0].set(x)
                                    
            endstep = time.time()
            
            # if info > 0:
            #     print("\nSTEP --", step, ": CONVERGED TO TOLERANCE NOT ACHIEVED")
            # elif info < 0:
            #     print("\nSTEP --", step, ": ILLEGAL INPUT OR BREAKDOWN")
            # else:
            #     print("\nSTEP --", step, ": SUCCESSFUL CONVERGED")
            
            print("\nSOLVE STEP " + str(step), "\ TIME SPEND: ", endstep - startstep, " SEC")
            
            U1[freedof, 0] += U0[freedof, 0]
            U[freedof, step] = U1[freedof, 0]
            U0[freedof, 0] = U1[freedof, 0]

            # U1.at[freedof, 0].add(U0[freedof, 0]) 
            # U.at[freedof, 0].set(U1[freedof, 0])  
            # U0.at[freedof, 0].set(U1[freedof, 0]) 
            
            # if "TRACKER" in solverset.keys():
            #     if solverset["TRACKER"]["show"]:
            #         plotset["step"] = step + 1
            #         plotset["val_list"] = U1
            #         plotset["fignumb"] = 1
            #         postprocset["TRACKER"] = solverset["TRACKER"]
            #         tracker_plot(
            #             postprocset, plotset, solverset["coord"], solverset["nodedof"]
            #         )
        
        solution['U'] = U
        return solution