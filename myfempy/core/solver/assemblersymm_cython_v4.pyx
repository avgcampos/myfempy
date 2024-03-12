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

@cython.cdivision(True)
@cython.exceptval(check=False)
@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function           
@cython.nonecheck(False)     
def getMatrixAssemblerSym(Model,
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

    mat = np.zeros((elemdof_max, elemdof_max), dtype=FLT64)
    cdef FLT64_t [:, ::1] mat_view = mat
    cdef FLT64_t VAL

    loc = np.zeros((elemdof_max,), dtype=INT32)
    cdef INT32_t [::1] loc_view = loc
    
    cdef INT32_t KI, KJ
    cdef INT32_t dim_band, dim_diag
    cdef Py_ssize_t ee, ii, jj, nb, nd

    dim_band = int(0.5*(elemdof_max*elemdof_max - elemdof_max)*elem_max)
    dim_diag = int(elemdof_max*elem_max)

    ith_band = np.zeros((dim_band,), dtype=INT32)
    cdef INT32_t [::1] ith_band_view = ith_band
    
    jth_band = np.zeros((dim_band,), dtype=INT32)
    cdef INT32_t [::1] jth_band_view = jth_band 
    
    val_band = np.zeros((dim_band,), dtype=FLT64)
    cdef FLT64_t [::1] val_band_view = val_band
    
    ith_diag = np.zeros((dim_diag,), dtype=INT32)
    cdef INT32_t [::1] ith_diag_view = ith_diag
    
    val_diag = np.zeros((dim_diag,), dtype=FLT64)
    cdef FLT64_t [::1] val_diag_view  = val_diag

    nb=0
    nd=0
    for ee in range(elem_max):
        mat_view = getMatrix(Model, inci, coord, tabmat, tabgeo, intgauss, ee, type_assembler)
        loc_view = getLoc(Model, inci, ee)
        for ii in range(elemdof_max):
            KI = loc_view[ii]
            for jj in range(ii, elemdof_max):
                KJ = loc_view[jj]
                VAL = mat_view[ii, jj]
                if KI == KJ:
                    ith_diag_view[nd]=KI
                    val_diag_view[nd]=VAL
                    nd+=1
                else:
                    ith_band_view[nb]=KI
                    jth_band_view[nb]=KJ
                    val_band_view[nb]=VAL
                    nb+=1

    return ith_diag, val_diag, ith_band, jth_band, val_band