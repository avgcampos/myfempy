cimport openmp
from cython cimport boundscheck, wraparound
from cython.parallel import parallel, prange

cimport numpy as np

ctypedef np.int32_t INT32_t
ctypedef np.float64_t FLT64_t

# cdef INT32_t index(INT32_t base_index, INT32_t elemdof, INT32_t ii, INT32_t jj) nogil:
#     idx = base_index + elemdof * ii + jj
#     return idx

@boundscheck(False) # turn off bounds-checking for entire function
@wraparound(False)  # turn off negative index wrapping for entire function  
def getVectorizationFull(INT32_t [::1] ith, INT32_t [::1] jth, FLT64_t [::1] val, INT32_t [::1] loc, FLT64_t [:, ::1] matrix, INT32_t ee, INT32_t elemdof, INT32_t num_threads):
    cdef Py_ssize_t ii, jj
    cdef Py_ssize_t LOOP_MAX = elemdof  
    cdef FLT64_t VAL
    cdef INT32_t KI, KJ
    cdef FLT64_t [:, ::1] mat_view = matrix
    cdef INT32_t [::1] loc_view = loc
    cdef INT32_t [::1] ith_view = ith
    cdef INT32_t [::1] jth_view = jth 
    cdef FLT64_t [::1] val_view = val
    cdef INT32_t idx
    cdef INT32_t size = elemdof * elemdof
    cdef INT32_t base_index = size * ee
    # with nogil, parallel(num_threads=num_threads):
    for ii in range(LOOP_MAX):
        KI = loc_view[ii]
        for jj in range(LOOP_MAX):
            KJ = loc_view[jj]
            VAL = mat_view[ii, jj]
            idx = base_index + elemdof * ii + jj
            ith_view[idx] = KI
            jth_view[idx] = KJ
            val_view[idx] = VAL
    return ith, jth, val
