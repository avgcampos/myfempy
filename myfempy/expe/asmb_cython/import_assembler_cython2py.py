import numpy as np
from scipy.sparse import coo_matrix

from myfempy.expe.asmb_cython.assembler_cython import getMatrixAssemblerSym_cy_v2

# @profile
def getMatrixAssemblerSYMM(Model, inci, coord, tabmat, tabgeo, intgauss, type_assembler):
        
    
    elem_set = Model.element.getElementSet()
    nodedof = len(elem_set["dofs"]['d'])

    shape_set = Model.shape.getShapeSet()
    nodecon = len(shape_set['nodes'])

    elemdof = nodecon * nodedof

    elemtot = inci.shape[0]
    nodetot = coord.shape[0]

    sdof = nodedof * nodetot
    
    rowsd, colsd, datad, rowsb, colsb, datab = getMatrixAssemblerSym_cy_v2(Model, inci, coord, tabmat, tabgeo, elemdof,  intgauss, type_assembler)
    
    # rowsd = data[0]
    # colsd = data[1]
    # datad = data[2]
    # rowsb = data[3]
    # colsb = data[4]
    # datab = data[5]
    
    mtKG_sp_sym = coo_matrix((datab, (rowsb, colsb)), shape=(sdof, sdof))
    mtKG_sp_sym += mtKG_sp_sym.transpose()
    mtKG_sp_sym += coo_matrix((datad, (rowsd, colsd)), shape=(sdof, sdof))
    
    return mtKG_sp_sym