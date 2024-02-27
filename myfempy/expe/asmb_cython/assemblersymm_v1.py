
import cython
from cython.parallel import prange


# import numpy as np
from scipy.sparse import coo_matrix, csc_matrix

# from myfempy.core.alglin import linsolve_spsolve
from myfempy.core.solver.assembler import Assembler, setAssembler
from myfempy.core.solver.assemblerfull import AssemblerFULL
# from myfempy.expe.asmb_cython.import_assembler_cython2py import getMatrixAssemblerSYMM

class AssemblerSYMMCY(Assembler):

    """
     Static Linear Solver Class <ConcreteClassService>
    """
    # @profile
    def getMatrixAssembler(Model, inci, coord, tabmat, tabgeo, intgauss, type_assembler):
        """
        getMatrixAssembler Assembler Symmertric Matrix System module <ConcreteClassService>

        Returns:
            matrix sparse
        """
        
        nodedof = cython.declare(cython.int)
        nodecon = cython.declare(cython.int)
        nodetot = cython.declare(cython.int)
        elemdof = cython.declare(cython.int)
        elemtot = cython.declare(cython.int)
        sdof = cython.declare(cython.int)
        ee = cython.declare(cython.int)
        ii = cython.declare(cython.int)
        jj = cython.declare(cython.int)
                
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
            mat, loc = setAssembler(Model, inci, coord, tabmat, tabgeo, intgauss, ee, type_assembler)
            for i in range(elemdof):
                for j in range(i, elemdof):
                    KI = loc[i]
                    KJ = loc[j]
                    val = mat[i, j]
                    if KI == KJ:
                        ith_diag.append(KI)
                        jth_diag.append(KJ)
                        val_diag.append(val)
                    else:
                        ith_band.append(KI)
                        jth_band.append(KJ)
                        val_band.append(val)
        A_sp_scipy = coo_matrix((val_band, (ith_band, jth_band)), shape=(sdof, sdof))
        A_sp_scipy += A_sp_scipy.transpose()
        A_sp_scipy += coo_matrix((val_diag, (ith_diag, jth_diag)), shape=(sdof, sdof))
        return A_sp_scipy
    
    def getLoadAssembler(loadaply, nodetot, nodedof):
        return AssemblerFULL.getLoadAssembler(loadaply, nodetot, nodedof)
    
    def getConstrains(constrains, nodetot, nodedof):
        return AssemblerFULL.getConstrains(constrains, nodetot, nodedof)