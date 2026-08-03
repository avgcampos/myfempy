from cython cimport boundscheck, cdivision, exceptval, nonecheck, wraparound
from libc.math cimport fabs

import numpy as np
cimport numpy as np

np.import_array()

ctypedef np.int32_t INT32
ctypedef np.float64_t FLT64

# ==============================================================================
# FUNÇÕES GEOMÉTRICAS E MATEMÁTICAS CDEF (INTERNAS)
# ==============================================================================

@cdivision(True)
@exceptval(check=False)
@boundscheck(False)
@wraparound(False)
@nonecheck(False)
cdef inline np.ndarray[FLT64, ndim=2, mode="c"] MATN(FLT64 r0):
    cdef np.ndarray[FLT64, ndim=2, mode="c"] N = np.zeros((4, 12), dtype=np.float64)
    N[0, 0]  = 0.5 * (1.0 - r0)
    N[0, 6]  = 0.5 * (1.0 + r0)
    N[1, 1]  = 0.25 * (r0 * r0 * r0 - 3.0 * r0 + 2.0)
    N[1, 5]  = 0.25 * (0.5 * r0 * r0 * r0 - 0.5 * r0 * r0 - 0.5 * r0 + 0.5)
    N[1, 7]  = 0.25 * (-r0 * r0 * r0 + 3.0 * r0 + 2.0)
    N[1, 11] = 0.25 * (0.5 * r0 * r0 * r0 + 0.5 * r0 * r0 - 0.5 * r0 - 0.5)  
    N[2, 2]  = 0.25 * (r0 * r0 * r0 - 3.0 * r0 + 2.0)
    N[2, 4]  = 0.25 * (0.5 * r0 * r0 * r0 - 0.5 * r0 * r0 - 0.5 * r0 + 0.5)
    N[2, 8]  = 0.25 * (-r0 * r0 * r0 + 3.0 * r0 + 2.0)
    N[2, 10] = 0.25 * (0.5 * r0 * r0 * r0 + 0.5 * r0 * r0 - 0.5 * r0 - 0.5)
    N[3, 3]  = 0.5 * (1.0 - r0)
    N[3, 9]  = 0.5 * (1.0 + r0)
    return N

@cdivision(True)
@exceptval(check=False)
@boundscheck(False)
@wraparound(False)
@nonecheck(False)
cdef inline np.ndarray[FLT64, ndim=2, mode="c"] MATDIFFN(FLT64 r0):
    cdef np.ndarray[FLT64, ndim=2, mode="c"] dN = np.zeros((1, 6), dtype=np.float64)
    dN[0, 0] = -0.5
    dN[0, 3] = 0.5
    return dN

@cdivision(True)
@exceptval(check=False)
@boundscheck(False)
@wraparound(False)
@nonecheck(False)
cdef inline np.ndarray[FLT64, ndim=2, mode="c"] MATDIFFDIFFN(FLT64 r0):
    cdef np.ndarray[FLT64, ndim=2, mode="c"] ddN = np.zeros((4, 12), dtype=np.float64)
    ddN[0, 0]  = -0.5
    ddN[0, 6]  = 0.5
    ddN[1, 1]  = 0.25 * (6.0 * r0)
    ddN[1, 5]  = 0.25 * (3.0 * r0 - 1.0)
    ddN[1, 7]  = 0.25 * (-6.0 * r0)
    ddN[1, 11] = 0.25 * (3.0 * r0 + 1.0)
    ddN[2, 2]  = 0.25 * (6.0 * r0)
    ddN[2, 4]  = 0.25 * (3.0 * r0 - 1.0)
    ddN[2, 8]  = 0.25 * (-6.0 * r0)
    ddN[2, 10] = 0.25 * (3.0 * r0 + 1.0)
    ddN[3, 3]  = -0.5
    ddN[3, 9]  = 0.5
    return ddN

@cdivision(True)
@exceptval(check=False)
@boundscheck(False)
@wraparound(False)
@nonecheck(False)
cdef inline FLT64 INV(FLT64 [:, ::1] A):
    return 1.0 / A[0, 0]

@cdivision(True)
@exceptval(check=False)
@boundscheck(False)
@wraparound(False)
@nonecheck(False)
cdef inline np.ndarray[FLT64, ndim=2, mode="c"] JACOBIANO(FLT64 r0, FLT64 [:, ::1] element_coord):  
    cdef np.ndarray[FLT64, ndim=2, mode="c"] diffN = MATDIFFN(r0)  
    cdef np.ndarray[FLT64, ndim=2, mode="c"] jac = np.zeros((1, 1), dtype=np.float64)
    jac[0, 0] = (diffN[0, 0] * element_coord[0, 0] + 
                 diffN[0, 1] * element_coord[1, 0] + 
                 diffN[0, 2] * element_coord[2, 0] + 
                 diffN[0, 3] * element_coord[3, 0] + 
                 diffN[0, 4] * element_coord[4, 0] + 
                 diffN[0, 5] * element_coord[5, 0])
    return jac

# ==============================================================================
# FUNÇÕES EXPORTADAS (PYTHON / DEF)
# ==============================================================================

@boundscheck(False)
@wraparound(False)
def ShapeFunctions(FLT64 [::1] r_coord, INT32 nodedof, FLT64 detJ):
    cdef np.ndarray[FLT64, ndim=2, mode="c"] mat_N = MATN(r_coord[0])
    mat_N[1, 5]  = 2.0 * detJ * mat_N[1, 5]
    mat_N[1, 11] = 2.0 * detJ * mat_N[1, 11]
    mat_N[2, 4]  = 2.0 * detJ * mat_N[2, 4]
    mat_N[2, 10] = 2.0 * detJ * mat_N[2, 10]
    return mat_N

@boundscheck(False)
@wraparound(False)
def DiffShapeFuntion(FLT64 [::1] r_coord, INT32 nodedof):
    return MATDIFFN(r_coord[0])

@boundscheck(False)
@wraparound(False)
def DiffDiffShapeFuntion(FLT64 [::1] r_coord, INT32 nodedof, FLT64 detJ):
    cdef np.ndarray[FLT64, ndim=2, mode="c"] mat_diff_N = MATDIFFDIFFN(r_coord[0])
    mat_diff_N[1, 5]  = 2.0 * detJ * mat_diff_N[1, 5]
    mat_diff_N[1, 11] = 2.0 * detJ * mat_diff_N[1, 11]
    mat_diff_N[2, 4]  = 2.0 * detJ * mat_diff_N[2, 4]
    mat_diff_N[2, 10] = 2.0 * detJ * mat_diff_N[2, 10]
    return mat_diff_N

@boundscheck(False)
@wraparound(False)
def Jacobian(FLT64 [::1] r_coord, FLT64 [:, ::1] element_coord):  
    return JACOBIANO(r_coord[0], element_coord)

@boundscheck(False)
@wraparound(False)
def invJacobi(FLT64 [::1] r_coord, FLT64 [:, ::1] element_coord, INT32 nodedof):
    cdef np.ndarray[FLT64, ndim=2, mode="c"] Jac = JACOBIANO(r_coord[0], element_coord)
    cdef FLT64 invJ = INV(Jac)
    cdef np.ndarray[FLT64, ndim=2, mode="c"] mat_invJ = np.zeros((4, 4), dtype=np.float64)
    mat_invJ[0, 0] = invJ
    mat_invJ[1, 1] = invJ * invJ
    mat_invJ[2, 2] = mat_invJ[1, 1]
    mat_invJ[3, 3] = mat_invJ[0, 0]
    return mat_invJ

@boundscheck(False)
@wraparound(False)
def detJacobi(FLT64 [::1] r_coord, FLT64 [:, ::1] element_coord):
    cdef np.ndarray[FLT64, ndim=2, mode="c"] Jac = JACOBIANO(r_coord[0], element_coord)
    return Jac[0, 0]

@boundscheck(False)
@wraparound(False)
def compute_B(INT32 [:, ::1] H, FLT64 [:, ::1] invJ, FLT64 [:, ::1] diffN):
    cdef INT32 h_rows = H.shape[0]
    cdef INT32 h_cols = H.shape[1]
    cdef INT32 invJ_cols = invJ.shape[1]
    cdef INT32 diffN_cols = diffN.shape[1]
    
    cdef np.ndarray[FLT64, ndim=2, mode="c"] T = np.zeros((h_rows, invJ_cols), dtype=np.float64)
    cdef np.ndarray[FLT64, ndim=2, mode="c"] B = np.zeros((h_rows, diffN_cols), dtype=np.float64)
    cdef Py_ssize_t i, j, k
    
    # 1. Multiplicação: T = H @ invJ
    for i in range(h_rows):
        for j in range(invJ_cols):
            for k in range(h_cols):
                T[i, j] += H[i, k] * invJ[k, j]
                
    # 2. Multiplicação: B = T @ diffN
    for i in range(h_rows):
        for j in range(diffN_cols):
            for k in range(invJ_cols):
                B[i, j] += T[i, k] * diffN[k, j]
                
    return B

@boundscheck(False)
@wraparound(False)
def StifLinear(
    FLT64 [::1] pt_gauss_nodes,
    FLT64 [::1] weight_factor,
    INT32 intgauss,
    FLT64 [:, ::1] element_coord,
    INT32 elemdof,
    INT32 nodedof,
    INT32 [:, ::1] H,
    FLT64 [:, ::1] C,
):
    cdef Py_ssize_t ip
    cdef INT32 r, c, k, l
    cdef INT32 c_rows = C.shape[0]
    cdef FLT64 detJ, scale, term, val_temp

    cdef np.ndarray[FLT64, ndim=1, mode="c"] current_pt_view = np.zeros(1, dtype=np.float64)
    cdef np.ndarray[FLT64, ndim=2, mode="c"] K_elem_mat = np.zeros((elemdof, elemdof), dtype=np.float64)
    cdef np.ndarray[FLT64, ndim=2, mode="c"] diffN, invJ, B

    for ip in range(intgauss):
        current_pt_view[0] = pt_gauss_nodes[ip]   

        detJ = detJacobi(current_pt_view, element_coord)
        diffN = DiffDiffShapeFuntion(current_pt_view, nodedof, detJ)
        invJ = invJacobi(current_pt_view, element_coord, nodedof)
        
        B = compute_B(H, invJ, diffN)
        scale = fabs(detJ) * weight_factor[ip]
        
        for k in range(c_rows):
            for l in range(c_rows):
                if C[k, l] != 0.0:
                    term = C[k, l] * scale
                    for r in range(elemdof):
                        if B[k, r] != 0.0:
                            val_temp = B[k, r] * term
                            for c in range(elemdof):
                                K_elem_mat[r, c] += val_temp * B[l, c]

    return K_elem_mat

@boundscheck(False)
@wraparound(False)
def MassLinear(
    FLT64 [::1] pt_gauss_nodes,
    FLT64 [::1] weight_factor,
    INT32 intgauss,
    FLT64 [:, ::1] element_coord,
    INT32 elemdof,
    INT32 nodedof,
    FLT64 [:, ::1] R,
):
    cdef Py_ssize_t ip
    cdef INT32 r, c, k, l
    cdef FLT64 detJ, scale, term, val_temp
    
    cdef np.ndarray[FLT64, ndim=1, mode="c"] current_pt_view = np.zeros(1, dtype=np.float64)
    cdef np.ndarray[FLT64, ndim=2, mode="c"] M_elem_mat = np.zeros((elemdof, elemdof), dtype=np.float64)
    cdef np.ndarray[FLT64, ndim=2, mode="c"] N

    for ip in range(intgauss):
        current_pt_view[0] = pt_gauss_nodes[ip]
        
        detJ = detJacobi(current_pt_view, element_coord)
        N = ShapeFunctions(current_pt_view, nodedof, detJ)
        
        scale = fabs(detJ) * weight_factor[ip]
        
        for k in range(4):
            for l in range(4):
                if R[k, l] != 0.0:
                    term = R[k, l] * scale
                    for r in range(elemdof):
                        if N[k, r] != 0.0:
                            val_temp = N[k, r] * term
                            for c in range(elemdof):
                                M_elem_mat[r, c] += val_temp * N[l, c]

    return M_elem_mat

@boundscheck(False)
@wraparound(False)
def NodeList(INT32 [:, ::1] inci, INT32 element_number):
    cdef np.ndarray[INT32, ndim=1, mode="c"] node_list = np.zeros(2, dtype=np.int32)
    node_list[0] = inci[element_number, 4]
    node_list[1] = inci[element_number, 5]
    return node_list
            
@boundscheck(False)
@wraparound(False)
def NodeCoord(FLT64 [:, ::1] coord, INT32 [::1] node_list):
    cdef np.ndarray[FLT64, ndim=2, mode="c"] element_coord = np.zeros((6, 1), dtype=np.float64)
    element_coord[0, 0] = coord[node_list[0] - 1, 1]
    element_coord[1, 0] = coord[node_list[0] - 1, 2]
    element_coord[2, 0] = coord[node_list[0] - 1, 3]
    element_coord[3, 0] = coord[node_list[1] - 1, 1]
    element_coord[4, 0] = coord[node_list[1] - 1, 2]
    element_coord[5, 0] = coord[node_list[1] - 1, 3]
    return element_coord

@boundscheck(False)
@wraparound(False)
def LocKey(INT32 [::1] node_list, INT32 nodedof):
    cdef np.ndarray[INT32, ndim=1, mode="c"] shape_key = np.zeros(2 * nodedof, dtype=np.int32)
    cdef Py_ssize_t node, dof
    for node in range(2):
        for dof in range(nodedof):
            shape_key[nodedof * node + dof] = nodedof * node_list[node] - (nodedof - dof)
    return shape_key