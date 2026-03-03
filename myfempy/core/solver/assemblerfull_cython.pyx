# distutils: language=c
# cython: language_level=3
# distutils: extra_compile_args=-fopenmp
# distutils: extra_link_args=-fopenmp
cimport openmp
from cython cimport boundscheck, wraparound

from cython.parallel import parallel, prange

cimport numpy as np

ctypedef np.int32_t INT32_t
ctypedef np.float64_t FLT64_t

@boundscheck(False) # turn off bounds-checking for entire function
@wraparound(False)  # turn off negative index wrapping for entire function           
def getVectorization(INT32_t [::1] ith, INT32_t [::1] jth, FLT64_t [::1] val, INT32_t [::1] loc, FLT64_t [:, ::1] matrix, INT32_t ee, INT32_t elemdof):
    cdef Py_ssize_t LOOP_MAX = elemdof  
    cdef FLT64_t VAL
    cdef INT32_t KI, KJ
    cdef Py_ssize_t ii, jj
    cdef FLT64_t [:, ::1] mat_view = matrix
    cdef INT32_t [::1] loc_view = loc
    cdef INT32_t [::1] ith_view = ith
    cdef INT32_t [::1] jth_view = jth 
    cdef FLT64_t [::1] val_view = val
    
    # openmp.omp_set_dynamic(0)
    # with nogil, parallel(num_threads=openmp.omp_get_num_threads()):
    for ii in range(LOOP_MAX):
        KI = loc_view[ii]
        for jj in range(LOOP_MAX):
            KJ = loc_view[jj]
            VAL = mat_view[ii, jj]
            ith_view[(LOOP_MAX*LOOP_MAX)*ee+LOOP_MAX*ii+jj]=KI
            jth_view[(LOOP_MAX*LOOP_MAX)*ee+LOOP_MAX*ii+jj]=KJ
            val_view[(LOOP_MAX*LOOP_MAX)*ee+LOOP_MAX*ii+jj]=VAL
    return ith, jth, val