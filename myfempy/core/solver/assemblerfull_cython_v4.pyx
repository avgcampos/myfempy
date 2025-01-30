# distutils: language=c
# cython: language_level=3
cimport numpy as np
from cython cimport boundscheck, wraparound

# from cython.view cimport array as cvarray
# from libc.stdint cimport int32_t

# cimport cython
# cimport openmp
# import numpy as np


# DTYPE = np.int64
# ctypedef np.int64_t DTYPE_t

# INT32 = np.int32
# FLT64 = np.float64

ctypedef np.int32_t INT32_t
ctypedef np.float64_t FLT64_t

# @cython.cdivision(True)
# @cython.exceptval(check=False)
@boundscheck(False) # turn off bounds-checking for entire function
@wraparound(False)  # turn off negative index wrapping for entire function           
# @cython.nonecheck(False)  
def getVectorizationFull(INT32_t [::1] ith, INT32_t [::1] jth, FLT64_t [::1] val, INT32_t [::1] loc, FLT64_t [:, ::1] matrix, INT32_t ee, INT32_t elemdof, INT32_t num_threads):

    # cdef unsigned long int num_threads
    cdef Py_ssize_t elemdof_max = elemdof  
    cdef FLT64_t VAL
    cdef INT32_t KI, KJ

    cdef Py_ssize_t ii, jj

    cdef FLT64_t [:, ::1] mat_view = matrix

    cdef INT32_t [::1] loc_view = loc
            
    cdef INT32_t [::1] ith_view = ith

    cdef INT32_t [::1] jth_view = jth 
    
    cdef FLT64_t [::1] val_view = val
                
    for ii in range(elemdof):
        KI = loc_view[ii]
        for jj in range(elemdof):
            KJ = loc_view[jj]
            VAL = mat_view[ii, jj]
            ith_view[(elemdof*elemdof)*ee+elemdof*ii+jj]=KI
            jth_view[(elemdof*elemdof)*ee+elemdof*ii+jj]=KJ
            val_view[(elemdof*elemdof)*ee+elemdof*ii+jj]=VAL
                    
    return ith_view, jth_view, val_view