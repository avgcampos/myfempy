from __future__ import annotations

from numpy import zeros, float64
from scipy.sparse import coo_matrix

from myfempy.core.solver.assembler import Assembler, setAssembler
from myfempy.core.solver.assemblerfull import AssemblerFULL
from myfempy.core.solver.assembler_cython_v2 import getMatrixAssemblerSym_cy
# from myfempy.core.solver.assembler_cython_v3 import getMatrixAssemblerSym_cy_parpool

class AssemblerSYMM(Assembler):

    """
     Assembler Symmetric Banded System Class <ConcreteClassService>
    """
    # @profile
    def getMatrixAssembler(Model, inci, coord, tabmat, tabgeo, intgauss, type_assembler):
        elem_set = Model.element.getElementSet()
        nodedof = len(elem_set["dofs"]['d'])
        shape_set = Model.shape.getShapeSet()
        nodecon = len(shape_set['nodes'])
        elemdof = nodecon * nodedof
        elemtot = inci.shape[0]
        nodetot = coord.shape[0]
        sdof = nodedof * nodetot
        
        rowsd, colsd, datad, rowsb, colsb, datab = getMatrixAssemblerSym_cy(Model, inci, coord, tabmat, tabgeo, elemdof,  intgauss, type_assembler)
        # rowsd, colsd, datad, rowsb, colsb, datab = getMatrixAssemblerSym_cy_parpool(Model, inci, coord, tabmat, tabgeo, elemdof,  intgauss, type_assembler)
            
        mtKG_sp_sym = coo_matrix((datab, (rowsb, colsb)), shape=(sdof, sdof))
        mtKG_sp_sym += mtKG_sp_sym.transpose()
        mtKG_sp_sym += coo_matrix((datad, (rowsd, colsd)), shape=(sdof, sdof))
        return mtKG_sp_sym
    
    def getLoadAssembler(loadaply, nodetot, nodedof):
        return AssemblerFULL.getLoadAssembler(loadaply, nodetot, nodedof)
    
    def getConstrains(constrains, nodetot, nodedof):
        return AssemblerFULL.getConstrains(constrains, nodetot, nodedof)