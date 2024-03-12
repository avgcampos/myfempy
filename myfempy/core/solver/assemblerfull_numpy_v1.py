# import cython
import numpy as np

from myfempy.core.solver.assembler import getMatrix, getLoc

def getMatrixAssembler_Full(Model,
                            inci,
                            coord,
                            tabmat,
                            tabgeo,
                            elemdof,
                            intgauss,
                            type_assembler):
    
    ith = np.zeros((inci.shape[0] * (elemdof * elemdof)), dtype=int)
    jth = np.zeros((inci.shape[0] * (elemdof * elemdof)), dtype=int)
    val = np.zeros((inci.shape[0] * (elemdof * elemdof)), dtype=float)
    for ee in range(inci.shape[0]):
        mat = getMatrix(Model, inci, coord, tabmat, tabgeo, intgauss, ee, type_assembler)
        loc = getLoc(Model, inci, ee)
        ith[(elemdof*elemdof)*ee:(elemdof*elemdof)*(ee+1)] = np.tile(loc.reshape(1, elemdof).T, (1, elemdof)).flatten("F")
        jth[(elemdof*elemdof)*ee:(elemdof*elemdof)*(ee+1)] = np.transpose(np.tile(loc.reshape(1, elemdof).T, (1, elemdof))).flatten("F")
        val[(elemdof*elemdof)*ee:(elemdof*elemdof)*(ee+1)] = mat.flatten("F")   
    return ith, jth, val