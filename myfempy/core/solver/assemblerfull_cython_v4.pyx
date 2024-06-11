cimport cython
cimport openmp

import numpy as np

cimport numpy as np

# DTYPE = np.int64
# ctypedef np.int64_t DTYPE_t

INT32 = np.int32
FLT64 = np.float64

ctypedef np.int32_t INT32_t
ctypedef np.float64_t FLT64_t

from myfempy.core.solver.assembler import getLoc, getMatrix


@cython.cdivision(True)
@cython.exceptval(check=False)
@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function           
@cython.nonecheck(False)  
def getMatrixAssembler_Full(Model,
                                 INT32_t [:, ::1] inci,
                                 FLT64_t [:, ::1] coord,
                                 FLT64_t [:, ::1] tabmat,
                                 FLT64_t [:, ::1] tabgeo,
                                 INT32_t elemdof,
                                 INT32_t intgauss,
                                 str type_assembler):

    # cdef unsigned long int num_threads
    cdef Py_ssize_t elemdof_max = elemdof
    cdef Py_ssize_t elem_max = inci.shape[0]    
    cdef FLT64_t VAL
    cdef INT32_t KI, KJ

    cdef Py_ssize_t ee, ii, jj

    mat = np.zeros((elemdof_max, elemdof_max), dtype=FLT64)
    cdef FLT64_t [:, ::1] mat_view = mat

    loc = np.zeros((elemdof_max,), dtype=INT32)
    cdef INT32_t [::1] loc_view = loc
            
    ith = np.zeros((elem_max * (elemdof_max * elemdof_max)), dtype=INT32)
    cdef INT32_t [::1] ith_view = ith
    
    jth = np.zeros((elem_max * (elemdof_max * elemdof_max)), dtype=INT32)
    cdef INT32_t [::1] jth_view = jth 
    
    val = np.zeros((elem_max * (elemdof_max * elemdof_max)), dtype=FLT64)
    cdef FLT64_t [::1] val_view = val
        
    for ee in range(elem_max):
        
        mat_view = getMatrix(Model, inci, coord, tabmat, tabgeo, intgauss, ee, type_assembler)
        loc_view = getLoc(Model, inci, ee)
        
        for ii in range(elemdof):
            KI = loc_view[ii]
            for jj in range(elemdof):
                KJ = loc_view[jj]
                VAL = mat_view[ii, jj]
                ith_view[(elemdof*elemdof)*ee+elemdof*ii+jj]=KI
                jth_view[(elemdof*elemdof)*ee+elemdof*ii+jj]=KJ
                val_view[(elemdof*elemdof)*ee+elemdof*ii+jj]=VAL
                    
    return ith_view, jth_view, val_view