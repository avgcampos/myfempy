# cimport openmp
from cython cimport boundscheck, cdivision, exceptval, nonecheck, wraparound

# from cython.parallel import parallel, prange

import numpy as np

cimport numpy as np

DTYPE = np.float64
ctypedef np.int32_t INT32
ctypedef np.float64_t FLT64

@cdivision(True)
@exceptval(check=False)
@boundscheck(False) # turn off bounds-checking for entire function
@wraparound(False)  # turn off negative index wrapping for entire function           
@nonecheck(False)   
cdef FLT64 [:, ::1] MATN(FLT64 [::1] r_coord):
    cdef FLT64 r0 = r_coord[0]
    N = np.zeros((4,12), dtype=np.float64)
    cdef FLT64 [:, ::1] N_view = N
    N_view[0][0] = 0.5*(1 - r0)
    N_view[0][6] = 0.5*(1 + r0)
    N_view[1][1] = 0.25*(r0**3 - 3*r0 + 2)
    N_view[1][5] = 0.25*(0.5*r0**3 - 0.5*r0**2 -0.5*r0 + 0.5)
    N_view[1][7] = 0.25*(r0**3 + 3*r0 + 2)
    N_view[1][11] = 0.25*(0.5*r0**3 + 0.5*r0**2 -0.5*r0 - 0.5)
    N_view[2][2] = 0.25*(r0**3 - 3*r0 + 2)
    N_view[2][4] = 0.25*(0.5*r0**3 - 0.5*r0**2 -0.5*r0 + 0.5)
    N_view[2][8] = 0.25*(-r0**3 + 3*r0 + 2)
    N_view[2][10] = 0.25*(0.5*r0**3 + 0.5*r0**2 -0.5*r0 - 0.5)
    N_view[3][3] = 0.5*(1 - r0)
    N_view[3][9] = 0.5*(1 + r0)
    return N

cdef FLT64 [:, ::1] MATDIFFN(FLT64 [::1] r):
    cdef FLT64 r0 = r[0]
    dN = np.zeros((1,6), dtype=np.float64)
    cdef FLT64 [:, ::1] dN_view = dN
    dN_view[0][0] = -0.5
    dN_view[0][3] = 0.5
    return dN

@cdivision(True)
@exceptval(check=False)
@boundscheck(False) # turn off bounds-checking for entire function
@wraparound(False)  # turn off negative index wrapping for entire function           
@nonecheck(False)   
cdef FLT64 [:, ::1] MATDIFFDIFFN(FLT64 [::1] r):
    cdef FLT64 r0 = r[0]
    ddN = np.zeros((4,12), dtype=np.float64)
    cdef FLT64 [:, ::1] ddN_view = ddN
    ddN_view[0][0] = -0.5
    ddN_view[0][6] = 0.5
    ddN_view[1][1] = 0.25*(6*r0)
    ddN_view[1][5] = 0.25*(3*r0 - 1)
    ddN_view[1][7] = 0.25*(-6*r0)
    ddN_view[1][11] = 0.25*(3*r0 + 1)
    ddN_view[2][2] = 0.25*(6*r0)
    ddN_view[2][4] = 0.25*(3*r0 - 1)
    ddN_view[2][8] = 0.25*(-6*r0)
    ddN_view[2][10] = 0.25*(3*r0 + 1)
    ddN_view[3][3] = -0.5
    ddN_view[3][9] = 0.5
    return ddN

# @cdivision(True)
# @exceptval(check=False)
# @boundscheck(False) # turn off bounds-checking for entire function
# @wraparound(False)  # turn off negative index wrapping for entire function           
# @nonecheck(False)   
# cdef FLT64 DET(FLT64 [:, ::1] A):
#     cdef FLT64 det = A[0][0]
#     return det

@cdivision(True)
@exceptval(check=False)
@boundscheck(False) # turn off bounds-checking for entire function
@wraparound(False)  # turn off negative index wrapping for entire function           
@nonecheck(False)   
cdef FLT64 INV(FLT64 [:, ::1] A):
    cdef FLT64 detA = A[0][0]
    cdef FLT64 invA = 1.0 / A[0][0]
    return invA

@cdivision(True)
@exceptval(check=False)
@boundscheck(False) # turn off bounds-checking for entire function
@wraparound(False)  # turn off negative index wrapping for entire function           
@nonecheck(False)   
cdef FLT64 [:, ::1] JACOBIANO(FLT64 [::1] r_coord, FLT64 [:, ::1] element_coord):  
    cdef FLT64 [:, ::1] diffN = MATDIFFN(r_coord)  
    jac = np.zeros((1, 1), dtype=np.float64)
    cdef FLT64 [:, ::1] jac_view = jac
    jac_view[0][0] = diffN[0][0]*element_coord[0][0]+diffN[0][1]*element_coord[1][0]+diffN[0][2]*element_coord[2][0]+diffN[0][3]*element_coord[3][0]+diffN[0][4]*element_coord[4][0]+diffN[0][5]*element_coord[5][0]
    return jac

@boundscheck(False) # turn off bounds-checking for entire function
@wraparound(False)  # turn off negative index wrapping for entire function                  
def ShapeFunctions(FLT64 [::1] r_coord, INT32 nodedof):
    N = MATN(r_coord)
    cdef FLT64 [:, ::1] shape_function = N
    matN = np.zeros((4, 2*nodedof), dtype=np.float64) 
    cdef FLT64 [:, ::1] mat_N = matN
    mat_N = shape_function
    
    return mat_N

@boundscheck(False) # turn off bounds-checking for entire function
@wraparound(False)  # turn off negative index wrapping for entire function           
def DiffShapeFuntion(FLT64 [::1] r_coord, INT32 nodedof):
    diffN = MATDIFFN(r_coord)
    cdef FLT64 [:, ::1] diff_shape_function = diffN
    matdiffN = np.zeros((1, 6), dtype=np.float64) 
    cdef FLT64 [:, ::1] mat_diff_N = matdiffN
    mat_diff_N = diff_shape_function
    
    return mat_diff_N

@boundscheck(False) # turn off bounds-checking for entire function
@wraparound(False)  # turn off negative index wrapping for entire function                
def DiffDiffShapeFuntion(FLT64 [::1] r_coord, INT32 nodedof):
    diffN = MATDIFFDIFFN(r_coord)
    cdef FLT64 [:, ::1] diff_shape_function = diffN
    matdiffN = np.zeros((4, 2*nodedof), dtype=np.float64) 
    cdef FLT64 [:, ::1] mat_diff_N = matdiffN
    mat_diff_N = diff_shape_function
    
    return mat_diff_N
    
@boundscheck(False) # turn off bounds-checking for entire function
@wraparound(False)  # turn off negative index wrapping for entire function           
def Jacobian(FLT64 [::1] r_coord, FLT64 [:, ::1] element_coord):  
    cdef FLT64 [:, ::1] Jac = JACOBIANO(r_coord, element_coord)
    return Jac
    
@boundscheck(False) # turn off bounds-checking for entire function
@wraparound(False)  # turn off negative index wrapping for entire function           
def invJacobi(FLT64 [::1] r_coord, FLT64 [:, ::1] element_coord, INT32 nodedof):
    cdef FLT64 [:, ::1] Jac = JACOBIANO(r_coord, element_coord)
    cdef FLT64 invJ = INV(Jac)
    mat_invJ = np.zeros((4, 4), dtype=np.float64)  
    cdef FLT64 [:, ::1] mat_invJ_view = mat_invJ
    mat_invJ_view[0][0] = invJ
    mat_invJ_view[1][1] = invJ*invJ
    mat_invJ_view[2][2] = mat_invJ_view[1][1]
    mat_invJ_view[3][3] = mat_invJ_view[0][0]
    return mat_invJ

@boundscheck(False) # turn off bounds-checking for entire function
@wraparound(False)  # turn off negative index wrapping for entire function           
def detJacobi(FLT64 [::1] r_coord, FLT64 [:, ::1] element_coord):
    cdef FLT64 [:, ::1] Jac = JACOBIANO(r_coord, element_coord)
    cdef FLT64 detJ = Jac[0][0]
    return detJ

@boundscheck(False) # turn off bounds-checking for entire function
@wraparound(False)  # turn off negative index wrapping for entire function           
def NodeList(INT32 [:, ::1] inci, INT32 element_number):
    cdef INT32 noi = inci[element_number, 4]
    cdef INT32 noj = inci[element_number, 5]
    cdef INT32 [::1] node_list = np.array([noi, noj], dtype=np.int32)                  
    return node_list
            
@boundscheck(False) # turn off bounds-checking for entire function
@wraparound(False)  # turn off negative index wrapping for entire function             
def NodeCoord(FLT64 [:, ::1] coord, INT32 [::1] node_list):
    cdef INT32 noi = node_list[0]
    cdef INT32 noj = node_list[1]
    cdef FLT64 xi = coord[noi - 1, 1]
    cdef FLT64 yi = coord[noi - 1, 2]
    cdef FLT64 zi = coord[noi - 1, 3]
    cdef FLT64 xj = coord[noj - 1, 1]
    cdef FLT64 yj = coord[noj - 1, 2]
    cdef FLT64 zj = coord[noj - 1, 3]
    cdef FLT64 [:, ::1] element_coord = np.array([[xi], [yi], [zi], [xj], [yj], [zj]], dtype=np.float64)
    return element_coord

@boundscheck(False) # turn off bounds-checking for entire function
@wraparound(False)  # turn off negative index wrapping for entire function           
def LocKey(INT32 [::1] node_list, INT32 nodedof):
    shape_key = np.zeros(2*nodedof, dtype=np.int32)
    cdef INT32 [::1] shape_key_view = shape_key
    cdef Py_ssize_t node, dof
    for node in range(len(node_list)):
        for dof in range(nodedof):
            shape_key_view[nodedof*node+dof] = nodedof * node_list[node] - (nodedof-dof)
    return shape_key