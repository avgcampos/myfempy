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
def getVectorization(INT32_t [::1] ith_band,
                           INT32_t [::1] jth_band,
                           FLT64_t [::1] val_band,
                           INT32_t [::1] ith_diag,
                           FLT64_t [::1] val_diag,
                           INT32_t nb,
                           INT32_t nd,
                           INT32_t [::1] loc,
                           FLT64_t [:, ::1] matrix,
                           INT32_t ee,
                           INT32_t elemdof):

    # cdef unsigned long int num_threads
    cdef Py_ssize_t LOOP_MAX = elemdof
    cdef Py_ssize_t ii, jj

    cdef FLT64_t [:, ::1] mat_view = matrix
    
    cdef INT32_t [::1] loc_view = loc
    
    cdef INT32_t KI, KJ
    cdef FLT64_t VAL
    

    cdef INT32_t [::1] ith_band_view = ith_band
    
    cdef INT32_t [::1] jth_band_view = jth_band 
    
    cdef FLT64_t [::1] val_band_view = val_band
    
    cdef INT32_t [::1] ith_diag_view = ith_diag
    
    cdef FLT64_t [::1] val_diag_view  = val_diag

    # cdef Py_ssize_t nb_view = nb
    # cdef Py_ssize_t nd_view = nd

    # openmp.omp_set_dynamic(0)
    # with nogil, parallel(num_threads=openmp.omp_get_num_threads()):
    for ii in range(LOOP_MAX):
        KI = loc_view[ii]
        for jj in range(ii, LOOP_MAX):
            KJ = loc_view[jj]
            VAL = mat_view[ii, jj]
            if KI == KJ:
                ith_diag_view[nd]=KI
                val_diag_view[nd]=VAL
                nd=nd+1
            else:
                ith_band_view[nb]=KI
                jth_band_view[nb]=KJ
                val_band_view[nb]=VAL
                nb=nb+1

    return ith_diag, val_diag, ith_band, jth_band, val_band, nb, nd