cimport cython
import numpy as np
cimport numpy as np
from cython.parallel cimport prange
# from cython.cimports.openmp import omp_set_dynamic, omp_get_num_threads

from myfempy.core.solver.assembler import setAssembler

@cython.cdivision(True)
@cython.exceptval(check=False)
@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
def getMatrixAssemblerSym_cy_parpool(Model,
                                 double[:, :] inci,
                                 double[:, :] coord,
                                 double[:, :] tabmat,
                                 double[:, :] tabgeo,
                                 int elemdof,
                                 int intgauss,
                                 str type_assembler):

    cdef unsigned long int i
    cdef unsigned long int j
    cdef unsigned long int num_threads
    cdef int element
    cdef list rowsb=[]
    cdef list colsb=[]
    cdef list datab=[]
    cdef list rowsd=[]
    cdef list colsd=[]
    cdef list datad=[]
    cdef list data=[]
    
    with nogil:
        for element in prange(inci.shape[0], schedule='dynamic'):
            with gil:
                data = setAssembler(Model, inci, coord, tabmat, tabgeo, intgauss, element, type_assembler)
            for i in prange(elemdof, schedule='dynamic'):
                for j in prange(i, elemdof, schedule='dynamic'):
                    with gil:
                        if data[1][i] == data[1][j]:
                            rowsd.append(data[1][i])
                            colsd.append(data[1][j])
                            datad.append(data[0][i, j])
                        else:
                            rowsb.append(data[1][i])
                            colsb.append(data[1][j])
                            datab.append(data[0][i, j])

    return np.array(rowsd, dtype=np.intp), np.array(colsd, dtype=np.intp), np.array(datad), \
           np.array(rowsb, dtype=np.intp), np.array(colsb, dtype=np.intp), np.array(datab)



# @cython.exceptval(check=False)
# @cython.boundscheck(False) # turn off bounds-checking for entire function
# @cython.wraparound(False)  # turn off negative index wrapping for entire function
# def getMatrixAssemblerSym_cy_v3(Model,
#                                 double[:, :] inci,
#                                 double[:, :] coord,
#                                 double[:, :] tabmat,
#                                 double[:, :] tabgeo,
#                                 int elemdof,
#                                 int intgauss):

#     cdef unsigned long int KI, KJ, i, j, elemtot, nodedof
#     cdef double [:, :] mat
#     cdef long [:] loc
#     cdef Py_ssize_t element
#     cdef list rowsb, colsb, datab, rowsd, colsd, datad
#     cdef list nodelist

#     elemtot = inci.shape[0]

#     rowsb = []
#     colsb = []
#     datab = []
#     rowsd = []
#     colsd = []
#     datad = []

#     elem_set = Model.element.getElementSet()
#     nodedof = len(elem_set["dofs"]['d'])

#     for element in prange(elemtot, nogil=True):
#         nodelist = Model.shape.getNodeList(inci, element)
#         mat = Model.element.getStifLinearMat(Model, inci, coord, tabmat, tabgeo, intgauss, element)
#         loc = Model.shape.getShapeKey(nodelist, nodedof)

#         for i in range(elemdof):
#             for j in range(i, elemdof):
#                 KI = loc[i]
#                 KJ = loc[j]
#                 val = mat[i, j]
#                 if KI == KJ:
#                     rowsd.append(KI)
#                     colsd.append(KJ)
#                     datad.append(val)
#                 else:
#                     rowsb.append(KI)
#                     colsb.append(KJ)
#                     datab.append(val)

#     return np.array(rowsd, dtype=np.intp), np.array(colsd, dtype=np.intp), np.array(datad), \
#            np.array(rowsb, dtype=np.intp), np.array(colsb, dtype=np.intp), np.array(datab)
