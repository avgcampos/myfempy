cimport cython
cimport openmp

from cython.parallel import parallel, prange

import numpy as np
cimport numpy as np

# DTYPE = np.int64
# ctypedef np.int64_t DTYPE_t

INT32 = np.int32
FLT64 = np.float64

ctypedef np.int32_t INT32_t
ctypedef np.float64_t FLT64_t

from myfempy.core.solver.assembler import getMatrix, getLoc


cdef FLT64_t [:,:] getMatrixCY(Model, INT32_t [:,:] inci, FLT64_t [:,:] coord, FLT64_t [:,:] tabmat, FLT64_t [:,:] tabgeo, INT32_t element, INT32_t intgauss, str type_assembler):
    return getMatrix(Model, inci, coord, tabmat, tabgeo, intgauss, element, type_assembler)


cdef INT32_t [:] getLocCY(Model, INT32_t [:,:] inci, INT32_t element):
    return getLoc(Model, inci, element)

cdef void getIndexVec(Model,
                     INT32_t [:,:] inci,
                     FLT64_t [:,:] coord,
                     FLT64_t [:,:] tabmat,
                     FLT64_t [:,:] tabgeo,
                     INT32_t elemdof, 
                     INT32_t intgauss, 
                     str type_assembler, 
                     INT32_t element, 
                     INT32_t [:] ith, 
                     INT32_t [:] jth,
                     FLT64_t [:] val):

    # # mat = np.zeros((elemdof_max, elemdof_max), dtype=FLT64)
    cdef FLT64_t [:, ::1] mat_view 

    # # loc = np.zeros((elemdof_max,), dtype=INT32)
    cdef INT32_t [::1] loc_view 

    with gil:
        mat_view = getMatrixCY(Model, inci, coord, tabmat, tabgeo, intgauss, element, type_assembler)
        loc_view = getLocCY(Model, inci, element)
    
    for ii in range(elemdof):
        KI = loc_view[ii]
        for jj in range(elemdof):
            KJ = loc_view[jj]
            VAL = mat_view[ii, jj]
            ith[(elemdof*elemdof)*element+elemdof*ii+jj]=KI
            jth[(elemdof*elemdof)*element+elemdof*ii+jj]=KJ
            val[(elemdof*elemdof)*element+elemdof*ii+jj]=VAL
    # return [ith, jth, val]



@cython.cdivision(True)
@cython.exceptval(check=False)
@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function           
@cython.nonecheck(False)  
def getMatrixAssemblerSym_cy_mp(Model,
                                 INT32_t [:,:] inci,
                                 FLT64_t [:,:] coord,
                                 FLT64_t [:,:] tabmat,
                                 FLT64_t [:,:] tabgeo,
                                 INT32_t elemdof,
                                 INT32_t intgauss,
                                 str type_assembler,
                                 INT32_t num_threads):

    # cdef unsigned long int num_threads
    cdef Py_ssize_t elemdof_max = elemdof
    cdef Py_ssize_t elem_max = inci.shape[0]    
    # cdef FLT64_t VAL
    # cdef INT32_t KI, KJ

    # cdef Py_ssize_t element, ii, jj

    # # mat = np.zeros((elemdof_max, elemdof_max), dtype=FLT64)
    # cdef FLT64_t [:, ::1] mat_view 

    # # loc = np.zeros((elemdof_max,), dtype=INT32)
    # cdef INT32_t [::1] loc_view 
            
    ith = np.zeros((elem_max * (elemdof_max * elemdof_max)), dtype=INT32)
    cdef INT32_t [::1] ith_view = ith
    
    jth = np.zeros((elem_max * (elemdof_max * elemdof_max)), dtype=INT32)
    cdef INT32_t [::1] jth_view = jth
    
    val = np.zeros((elem_max * (elemdof_max * elemdof_max)), dtype=FLT64)
    cdef FLT64_t [::1] val_view = val

    # assert tuple(mat_view) == tuple(loc_view)

    with nogil, parallel(num_threads=num_threads):
        
        for element in prange(elem_max, schedule='static'):
            
            getIndexVec(Model, inci, coord, tabmat, tabgeo, elemdof_max, intgauss, type_assembler, element, ith_view, jth_view, val_view)

            # with gil:
            #     mat_view = getMatrixCY(Model, inci, coord, tabmat, tabgeo, intgauss, ee, type_assembler)
            #     loc_view = getLocCY(Model, inci, ee)
            
            # for ii in range(elemdof_max):
            #     KI = loc_view[ii]
            #     for jj in range(elemdof_max):
            #         KJ = loc_view[jj]
            #         VAL = mat_view[ii, jj]
            #         ith_view[(elemdof_max*elemdof_max)*ee+elemdof_max*ii+jj]=KI
            #         jth_view[(elemdof_max*elemdof_max)*ee+elemdof_max*ii+jj]=KJ
            #         val_view[(elemdof_max*elemdof_max)*ee+elemdof_max*ii+jj]=VAL
                    
    return ith, jth, val