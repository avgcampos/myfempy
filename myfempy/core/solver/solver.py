from abc import ABC, abstractmethod
from re import I

import numpy as np
import scipy.sparse as sp

# def getSolver(set_solver):
#     if set_solver['solver'] == 'SLI':
#         from myfempy.core.staticlinear import StaticLinear
#         return StaticLinear
#     elif set_solver['solver'] == 'EIG':
#         from myfempy.core.modallinear import ModalLinear
#         return ModalLinear
#     else:
#         pass

class Solver(ABC):
    
    """
     Solver API Class <ClassService>
    """
        
    @abstractmethod
    def runSolve():
        
        """
        Solve Solve FE System: Model + Physics
        """
        pass
    
    
    def getMatrixAssemblerFULL(Model, inci, coord, tabmat, tabgeo, intgauss, type_assembler):
        
        """
        getMatrixAssembler Assembler module <ConcreteClassService>

        Returns:
            matrix sparse
        """    
        elem_set = Model.element.getElementSet()
        nodedof = len(elem_set["dofs"]['d'])
        shape_set = Model.shape.getShapeSet()
        nodecon = len(shape_set['nodes'])
        elemdof = nodecon * nodedof
        elemtot = len(inci)
        nodetot = len(coord)
        ith = np.zeros((elemtot * (elemdof * elemdof)), dtype=int)
        jth = np.zeros((elemtot * (elemdof * elemdof)), dtype=int)
        val = np.zeros((elemtot * (elemdof * elemdof)), dtype=float)
        for ee in range(elemtot):
            mat, loc = Solver.__setassembler(Model, inci, coord, tabmat, tabgeo, intgauss, ee, type_assembler)
            ith[(elemdof*elemdof)*ee:(elemdof*elemdof)*(ee+1)] = np.tile(loc.reshape(1, elemdof).T, (1, elemdof)).flatten("F")
            jth[(elemdof*elemdof)*ee:(elemdof*elemdof)*(ee+1)] = np.transpose(np.tile(loc.reshape(1, elemdof).T, (1, elemdof))).flatten("F")
            val[(elemdof*elemdof)*ee:(elemdof*elemdof)*(ee+1)] = mat.flatten("F")       
        A_sp_scipy_csc = sp.csc_matrix((val, (ith, jth)), shape=(nodedof * nodetot, nodedof * nodetot),)        
        return A_sp_scipy_csc

    
    def getMatrixAssemblerSYMM(Model, inci, coord, tabmat, tabgeo, intgauss, type_assembler):
        
        """
        getMatrixAssembler Assembler Symmertric Matrix System module <ConcreteClassService>

        Returns:
            matrix sparse
        """
                
        elem_set = Model.element.getElementSet()
        nodedof = len(elem_set["dofs"]['d'])

        shape_set = Model.shape.getShapeSet()
        nodecon = len(shape_set['nodes'])

        elemdof = nodecon * nodedof

        elemtot = len(inci)
        nodetot = len(coord)

        sdof = nodedof * nodetot

        ith_band = []
        jth_band = []
        val_band = []
        ith_diag = []
        jth_diag = []
        val_diag = []
        
        for ee in range(elemtot):
            mat, loc = Solver.__setassembler(Model, inci, coord, tabmat, tabgeo, intgauss, ee, type_assembler)
            for ii in range(elemdof):
                KI = loc[ii]         
                for jj in range(elemdof):
                    KJ = loc[jj]
                    if KI==KJ:
                        ith_diag.append(KI)
                        jth_diag.append(KJ)
                        val_diag.append(mat[ii, jj])
                    elif KI<KJ:
                        ith_band.append(KI)
                        jth_band.append(KJ)
                        val_band.append(mat[ii, jj])
        A_sp_scipy = sp.coo_matrix((val_band, (ith_band, jth_band)), shape=(sdof, sdof))
        A_sp_scipy += A_sp_scipy.transpose()
        A_sp_scipy += sp.coo_matrix((val_diag, (ith_diag, jth_diag)), shape=(sdof, sdof))
        return A_sp_scipy
    
    def getLoadAssembler(loadaply, nodetot, nodedof):
        
        """
        getLoadAssembler Assembler module <ConcreteClassService>

        Returns:
            vector sparse
        """
        
        forcevec = np.zeros((nodedof * nodetot,len(np.unique(loadaply[:, 3])),))
        
        for fstep in range(len(np.unique(loadaply[:, 3]))):
            forceaply = loadaply[np.where(loadaply[:, 3] == fstep + 1), :][0]
            nload = forceaply.shape[0]
            
            if nodedof == 1:
                for ii in range(nload):
                    if int(forceaply[ii, 1]) == 1:
                        gdlload = int(nodedof * forceaply[ii, 0] - (nodedof- 1))
                        forcevec[gdlload, fstep] += forceaply[ii, 2]
            
            elif nodedof == 2:
                for ii in range(nload):
                    if int(forceaply[ii, 1]) == 1:
                        gdlload = int(nodedof * forceaply[ii, 0] - (nodedof))
                        forcevec[gdlload, fstep] += forceaply[ii, 2]
                    
                    elif int(forceaply[ii, 1]) == 2:
                        gdlload = int(nodedof * forceaply[ii, 0] - (nodedof - 1))
                        forcevec[gdlload, fstep] += forceaply[ii, 2]
                        
            elif nodedof == 3:
                for ii in range(nload):
                    if int(forceaply[ii, 1]) == 1:
                        gdlload = int(nodedof * forceaply[ii, 0] - (nodedof))
                        forcevec[gdlload, fstep] += forceaply[ii, 2]
                    
                    elif int(forceaply[ii, 1]) == 2:
                        gdlload = int(nodedof * forceaply[ii, 0] - (nodedof - 1))
                        forcevec[gdlload, fstep] += forceaply[ii, 2]
                        
                    elif int(forceaply[ii, 1]) == 3:
                        gdlload = int(nodedof * forceaply[ii, 0] - (nodedof - 2))
                        forcevec[gdlload, fstep] += forceaply[ii, 2]
        
            else:
                pass
        forcevec = sp.csc_matrix(forcevec)
        return forcevec
    

    def getConstrains(constrains, nodetot, nodedof):
        
        """
        getConstrains Constrain module <ConcreteClassService>

        Returns:
            _description_
        """
              
        ntbc = len(constrains)
        fixedof = np.zeros((1, nodedof * nodetot))
        
        if nodedof == 1:
            for ii in range(ntbc):
                no = int(constrains[ii, 1])
                if int(constrains[ii, 0]) == 1:
                    fixedof[0, nodedof * no - 1] = (nodedof * no)
        
        elif nodedof == 2:
            for ii in range(ntbc):
                no = int(constrains[ii, 1])
                if int(constrains[ii, 0]) == 0:
                    fixedof[0, nodedof * no - 2] = (nodedof * no - 1)
                    fixedof[0, nodedof * no - 1] = (nodedof * no)
                elif int(constrains[ii, 0]) == 1:
                    fixedof[0, nodedof * no - 2] = (nodedof * no - 1)
                elif int(constrains[ii, 0]) == 2:
                    fixedof[0, nodedof * no - 1] = (nodedof * no)
            
        elif nodedof == 3:
            for ii in range(ntbc):
                no = int(constrains[ii, 1])
                if int(constrains[ii, 0]) == 0:
                    fixedof[0, nodedof * no - 3] = (nodedof * no - 2)
                    fixedof[0, nodedof * no - 2] = (nodedof * no - 1)
                    fixedof[0, nodedof * no - 1] = (nodedof * no)
                elif int(constrains[ii, 0]) == 1:
                    fixedof[0, nodedof * no - 3] = (nodedof * no - 2)
                elif int(constrains[ii, 0]) == 2:
                    fixedof[0, nodedof * no - 2] = (nodedof * no - 1)
                elif int(constrains[ii, 0]) == 3:
                    fixedof[0, nodedof * no - 1] = (nodedof * no)
            
        fixedof = fixedof[np.nonzero(fixedof)]
        fixedof = fixedof - np.ones_like(fixedof)
        alldof = np.arange(0, nodedof * nodetot, 1, int)
        freedof = np.setdiff1d(alldof, fixedof)
        return freedof, fixedof
    
    def addMatrix(A, A_add, loc):
        """
        addMatrix sum/add two matrix

        Arguments:
            A -- matrix numpy array
            A_add -- matrix numpy array to sum with A
            loc -- dof relative

        Returns:
            matrix numpy array
        """
        
        A[np.ix_(loc, loc)] += A_add
        return A
    
    def setSteps(steps):
        
        """
        setSteps steps setting

        Arguments:
            steps -- dict

        Returns:
            nsteps -- Number of steps [int]
        """
        
        start = steps["start"]
        end = steps["end"]
        substep = steps["step"]
        if (end - start) == 0:
            nsteps = int(end)
        else:
            nsteps = int((end - start) / substep)
        return nsteps
    
        
    def getLinSysSolve(A, b):
        """
        linsolve_direct Solve the sparse linear system Ax=b, where b may be a
        vector or a matrix

        Arguments:
            A -- The square matrix A will be converted into CSC or CSR form
            b -- The matrix or vector representing the right hand side of the
            equation

        Returns:
            x -- The solution of the sparse linear equation
        """
        x = sp.linalg.spsolve(A, b)
        return x

    def getGenMinResSolve(A, b, x0, M):
        """
        linsolve_gmres Use Generalized Minimal RESidual iteration to solve Ax=b

        Arguments:
            A -- The real or complex N-by-N matrix of the linear system
            b -- Right hand side of the linear system
            x0 -- Starting guess for the solution (a vector of zeros by default)
            M -- Inverse of the preconditioner of A

        Returns:
            x -- The solution of the sparse linear equation
            info -- Provides convergence information:
                * 0  : successful exit
                * >0 : convergence to tolerance not achieved, number of iterations
                * <0 : illegal input or breakdown
        """
        x, info =  sp.linalg.gmres(A=A, b=b, x0=x0, tol=1E-10, M=M, maxiter=1000)
        return x, info

    def getEigHerSysSolve(A, M, k):
        """
        eigsolve_eigsh Find k eigenvalues and eigenvectors of the real symmetric
        square matrix or complex Hermitian matrix A to solves A@x[i]=w[i]*M@x[i]

        Arguments:
            A -- A square operator representing the operation A @ x, where A
            is real symmetric or complex Hermitian
            M -- The operation M @ x for the generalized eigenvalue problem
            A@x[i]=w[i]*M@x[i]
                For best results, the data type of M should be the same as
                that of A. Additionally:
                    If sigma is None, M is symmetric positive definite
                    If sigma is specified, M is symmetric positive semi-definite
                    In buckling mode, M is symmetric indefinite
            k -- The number of eigenvalues and eigenvectors desired

        Returns:
            w -- Array of k eigenvalues.
            v -- An array representing the k eigenvectors
        """
        w, v = sp.linalg.eigsh(A, k, M, sigma=1, which="LM", maxiter=1000)
        return w, v
        
        
    def __setassembler(Model, inci, coord, tabmat, tabgeo, intgauss, element_number, type_assembler):

        elem_set = Model.element.getElementSet()
        nodedof = len(elem_set["dofs"]['d'])
        
        nodelist = Model.shape.getNodeList(inci, element_number)
         
        if type_assembler == 'linear_stiffness':
            mat = Model.element.getStifLinearMat(Model, inci, coord, tabmat, tabgeo, intgauss, element_number)
            loc = Model.shape.getShapeKey(nodelist, nodedof)
            return mat, loc
        
        elif type_assembler == 'mass_consistent':
            mat = Model.element.getMassConsistentMat(Model, inci, coord, tabmat, tabgeo, intgauss, element_number)
            loc = Model.shape.getShapeKey(nodelist, nodedof)
            return mat, loc
        
        elif type_assembler == 'mass_lumped':
            mat = Model.element.getMassLumpedMat(Model, inci, coord, tabmat, tabgeo, intgauss, element_number)
            loc = Model.shape.getShapeKey(nodelist, nodedof)
            return mat, loc