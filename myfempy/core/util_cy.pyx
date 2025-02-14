import numpy as np
cimport numpy as np
from cython cimport boundscheck, wraparound, cdivision
from cython.parallel import prange, parallel

ctypedef np.float64_t REAL_t
ctypedef np.int32_t INT32
from scipy.linalg.blas import ddot, sdot  #, sdotd

"""
R = A * B (warning: A, B should't have been numpy.array.transposed!)
@R: result matrix
@n_jobs: num of CPUs to be used
"""
@boundscheck(False)  # Turn off bounds-checking for entire function
@wraparound(False)   # Turn off negative index wrapping for entire function
@cdivision(True)
def fast_dot(np.ndarray[REAL_t, ndim=2] A, np.ndarray[REAL_t, ndim=2] B, np.ndarray[REAL_t, ndim=2] R, int n_jobs):
    cdef INT32 A_row = <INT32>A.shape[0]
    cdef INT32 A_col = <INT32>A.shape[1]
    cdef INT32 B_row = <INT32>B.shape[0]
    cdef INT32 B_col = <INT32>B.shape[1]
    
    # Declare memoryviews
    cdef REAL_t *mat_1 = <REAL_t *>np.PyArray_DATA(A) #A.data
    cdef REAL_t *mat_2 = <REAL_t *>np.PyArray_DATA(B)
    cdef REAL_t *result = <REAL_t *>np.PyArray_DATA(R)
    cdef INT32 cores = <INT32>n_jobs

    cdef Py_ssize_t i, j, k

    # Use OpenMP parallelism with static scheduling for better load balancing
    with nogil, parallel(num_threads=cores):
        for i in prange(A_row, schedule='static'):
            for j in range(B_col):
                result[i*B_col+j] = <REAL_t>0.0
                for k in range(A_col):
                    result[i*B_col+j] += mat_1[i*A_col+k] * mat_2[k*B_col+j]


# @boundscheck(False) # turn off bounds-checking for entire function
# @wraparound(False)  # turn off negative index wrapping for entire function   
# @cdivision(True)
# def fast_dot(np.ndarray[REAL_t, ndim=2] A, np.ndarray[REAL_t, ndim=2] B, np.ndarray[REAL_t, ndim=2] R, int n_jobs):
#     cdef INT32 A_row = <INT32>A.shape[0]
#     cdef INT32 A_col = <INT32>A.shape[1]
#     cdef INT32 B_row = <INT32>B.shape[0]
#     cdef INT32 B_col = <INT32>B.shape[1]
#     cdef REAL_t *mat_1 = <REAL_t *>np.PyArray_DATA(A) #A.data
#     cdef REAL_t *mat_2 = <REAL_t *>np.PyArray_DATA(B)
#     cdef REAL_t *result = <REAL_t *>np.PyArray_DATA(R)
#     cdef INT32 cores = <INT32>n_jobs

#     cdef Py_ssize_t i, j, k
#     for i in prange(A_row, nogil=True, num_threads=cores, schedule='static'):
#         for j in range(B_col):
#             result[i*B_col+j] = <REAL_t>0.0
#             for k in range(A_col):
#                 result[i*B_col+j] += (&mat_1[i*A_col+k]) * (&mat_2[k*B_col+j])











# """
# R = A * B^T (warning: A, B should't have been numpy.array.transposed!)
# @R: result matrix
# @n_jobs: num of CPUs to be used
# """
# @boundscheck(False)
# @wraparound(False)
# @cdivision(True)
# def fast_dot_blas(np.ndarray[REAL_t, ndim=2] A, np.ndarray[REAL_t, ndim=2] B, np.ndarray[REAL_t, ndim=2] R, INT32 n_jobs):
#     cdef INT32 A_row = <int>A.shape[0]
#     cdef INT32 A_col = <int>A.shape[1]
#     cdef INT32 B_row = <int>B.shape[0]
#     cdef INT32 B_col = <int>B.shape[1]
#     cdef REAL_t *mat_1 = <REAL_t *>np.PyArray_DATA(A) #A.data
#     cdef REAL_t *mat_2 = <REAL_t *>np.PyArray_DATA(B)
#     cdef REAL_t *result = <REAL_t *>np.PyArray_DATA(R)
#     cdef int cores = <int>n_jobs

#     cdef INT32 ONE = 1
#     cdef Py_ssize_t i, j
#     for i in range(A_row):
#     # for i in prange(A_row, nogil=True, num_threads=cores, schedule='static'):
#         for j in range(B_row):
#             result[i*B_row+j] = <REAL_t>sdot(A_col, mat_1[i*A_col], ONE, mat_2[j*B_col], ONE)