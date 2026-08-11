from cython cimport boundscheck, cdivision, exceptval, nonecheck, wraparound
from libc.math cimport fabs

import numpy as np
cimport numpy as np

np.import_array()

ctypedef np.int32_t INT32
ctypedef np.float64_t FLT64

# ==============================================================================
# FUNÇÕES GEOMÉTRICAS E MATEMÁTICAS CDEF (INTERNAS - ARRAYS C PURO)
# ==============================================================================
    
@cdivision(True)
@exceptval(check=False)
@boundscheck(False)
@wraparound(False)         
@nonecheck(False)
cdef inline void MATN(FLT64 r0, FLT64 r1, FLT64 r2, FLT64 N[8]) nogil:
    N[0] = 0.125 * (1.0 - r0) * (1.0 - r1) * (1.0 - r2)
    N[1] = 0.125 * (1.0 + r0) * (1.0 - r1) * (1.0 - r2)
    N[2] = 0.125 * (1.0 + r0) * (1.0 + r1) * (1.0 - r2)
    N[3] = 0.125 * (1.0 - r0) * (1.0 + r1) * (1.0 - r2)
    N[4] = 0.125 * (1.0 - r0) * (1.0 - r1) * (1.0 + r2)
    N[5] = 0.125 * (1.0 + r0) * (1.0 - r1) * (1.0 + r2)
    N[6] = 0.125 * (1.0 + r0) * (1.0 + r1) * (1.0 + r2)
    N[7] = 0.125 * (1.0 - r0) * (1.0 + r1) * (1.0 + r2)
    
@cdivision(True)
@exceptval(check=False)
@boundscheck(False)
@wraparound(False)         
@nonecheck(False)
cdef inline void MATDIFFN(FLT64 r0, FLT64 r1, FLT64 r2, FLT64 dN[3][8]) nogil:
    dN[0][0] = 0.125 * (-1.0 + r1) * (1.0 - r2)
    dN[0][1] = 0.125 * (1.0 - r1) * (1.0 - r2)
    dN[0][2] = 0.125 * (1.0 + r1) * (1.0 - r2)
    dN[0][3] = 0.125 * (-1.0 - r1) * (1.0 - r2)
    dN[0][4] = 0.125 * (-1.0 + r1) * (1.0 + r2)
    dN[0][5] = 0.125 * (1.0 - r1) * (1.0 + r2)
    dN[0][6] = 0.125 * (1.0 + r1) * (1.0 + r2)
    dN[0][7] = 0.125 * (-1.0 - r1) * (1.0 + r2)
    dN[1][0] = 0.125 * (-1.0 + r0) * (1.0 - r2)
    dN[1][1] = 0.125 * (-1.0 - r0) * (1.0 - r2)
    dN[1][2] = 0.125 * (1.0 + r0) * (1.0 - r2)
    dN[1][3] = 0.125 * (1.0 - r0) * (1.0 - r2)
    dN[1][4] = 0.125 * (-1.0 + r0) * (1.0 + r2)
    dN[1][5] = 0.125 * (-1.0 - r0) * (1.0 + r2)
    dN[1][6] = 0.125 * (1.0 + r0) * (1.0 + r2)
    dN[1][7] = 0.125 * (1.0 - r0) * (1.0 + r2)
    dN[2][0] = 0.125 * (-1.0 + r0) * (1.0 - r1)
    dN[2][1] = 0.125 * (-1.0 - r0) * (1.0 - r1)
    dN[2][2] = 0.125 * (-1.0 - r0) * (1.0 + r1)
    dN[2][3] = 0.125 * (-1.0 + r0) * (1.0 + r1)
    dN[2][4] = 0.125 * (1.0 - r0) * (1.0 - r1)
    dN[2][5] = 0.125 * (1.0 + r0) * (1.0 - r1)
    dN[2][6] = 0.125 * (1.0 + r0) * (1.0 + r1)
    dN[2][7] = 0.125 * (1.0 - r0) * (1.0 + r1)
    
@cdivision(True)
@exceptval(check=False)
@boundscheck(False)
@wraparound(False)         
@nonecheck(False)
cdef inline void INV(FLT64 A[3][3], FLT64 invA[3][3]) nogil:
    cdef FLT64 detA = A[0][0] * A[1][1] * A[2][2] + A[0][1] * A[1][2] * A[2][0] + A[0][2] * A[1][0] * A[2][1] - A[2][0] * A[1][1] * A[0][2] - A[2][1] * A[1][2] * A[0][0] - A[2][2] * A[1][0] * A[0][1]
    cdef FLT64 invDet = 1.0 / detA
    invA[0][0] = invDet * (A[1][1] * A[2][2] - A[2][1] * A[1][2])
    invA[0][1] = -invDet * (A[0][1] * A[2][2] - A[2][1] * A[0][2])
    invA[0][2] = invDet * (A[0][1] * A[1][2] - A[1][1] * A[0][2])
    invA[1][0] = -invDet * (A[1][0] * A[2][2] - A[2][0] * A[1][2])
    invA[1][1] = invDet * (A[0][0] * A[2][2] - A[2][0] * A[0][2])
    invA[1][2] = -invDet * (A[0][0] * A[1][2] - A[1][0] * A[0][2])
    invA[2][0] = invDet * (A[1][0] * A[2][1] - A[2][0] * A[1][1])
    invA[2][1] = -invDet * (A[0][0] * A[2][1] - A[2][0] * A[0][1])
    invA[2][2] = invDet * (A[0][0] * A[1][1] - A[1][0] * A[0][1])

    
@cdivision(True)
@exceptval(check=False)
@boundscheck(False)
@wraparound(False)         
@nonecheck(False)
cdef inline void JACOBIANO(FLT64 r0, FLT64 r1, FLT64 r2, FLT64 [:, ::1] element_coord, FLT64 jac[3][3]) nogil:
    cdef FLT64 diffN[3][8]
    MATDIFFN(r0, r1, r2, diffN)
    
    jac[0][0] = (
        diffN[0][0]*element_coord[0, 0] + diffN[0][1]*element_coord[1, 0] +
        diffN[0][2]*element_coord[2, 0] + diffN[0][3]*element_coord[3, 0] +
        diffN[0][4]*element_coord[4, 0] + diffN[0][5]*element_coord[5, 0] +
        diffN[0][6]*element_coord[6, 0] + diffN[0][7]*element_coord[7, 0]
    )
    jac[0][1] = (
        diffN[0][0]*element_coord[0, 1] + diffN[0][1]*element_coord[1, 1] +
        diffN[0][2]*element_coord[2, 1] + diffN[0][3]*element_coord[3, 1] +
        diffN[0][4]*element_coord[4, 1] + diffN[0][5]*element_coord[5, 1] +
        diffN[0][6]*element_coord[6, 1] + diffN[0][7]*element_coord[7, 1]
    )
    jac[0][2] = (
        diffN[0][0]*element_coord[0, 2] + diffN[0][1]*element_coord[1, 2] +
        diffN[0][2]*element_coord[2, 2] + diffN[0][3]*element_coord[3, 2] +
        diffN[0][4]*element_coord[4, 2] + diffN[0][5]*element_coord[5, 2] +
        diffN[0][6]*element_coord[6, 2] + diffN[0][7]*element_coord[7, 2]
    )
    jac[1][0] = (
        diffN[1][0]*element_coord[0, 0] + diffN[1][1]*element_coord[1, 0] +
        diffN[1][2]*element_coord[2, 0] + diffN[1][3]*element_coord[3, 0] +
        diffN[1][4]*element_coord[4, 0] + diffN[1][5]*element_coord[5, 0] +
        diffN[1][6]*element_coord[6, 0] + diffN[1][7]*element_coord[7, 0]
    )
    jac[1][1] = (
        diffN[1][0]*element_coord[0, 1] + diffN[1][1]*element_coord[1, 1] +
        diffN[1][2]*element_coord[2, 1] + diffN[1][3]*element_coord[3, 1] +
        diffN[1][4]*element_coord[4, 1] + diffN[1][5]*element_coord[5, 1] +
        diffN[1][6]*element_coord[6, 1] + diffN[1][7]*element_coord[7, 1]
    )
    jac[1][2] = (
        diffN[1][0]*element_coord[0, 2] + diffN[1][1]*element_coord[1, 2] +
        diffN[1][2]*element_coord[2, 2] + diffN[1][3]*element_coord[3, 2] +
        diffN[1][4]*element_coord[4, 2] + diffN[1][5]*element_coord[5, 2] +
        diffN[1][6]*element_coord[6, 2] + diffN[1][7]*element_coord[7, 2]
    )
    jac[2][0] = (
        diffN[2][0]*element_coord[0, 0] + diffN[2][1]*element_coord[1, 0] +
        diffN[2][2]*element_coord[2, 0] + diffN[2][3]*element_coord[3, 0] +
        diffN[2][4]*element_coord[4, 0] + diffN[2][5]*element_coord[5, 0] +
        diffN[2][6]*element_coord[6, 0] + diffN[2][7]*element_coord[7, 0]
    )
    jac[2][1] = (
        diffN[2][0]*element_coord[0, 1] + diffN[2][1]*element_coord[1, 1] +
        diffN[2][2]*element_coord[2, 1] + diffN[2][3]*element_coord[3, 1] +
        diffN[2][4]*element_coord[4, 1] + diffN[2][5]*element_coord[5, 1] +
        diffN[2][6]*element_coord[6, 1] + diffN[2][7]*element_coord[7, 1]
    )
    jac[2][2] = (
        diffN[2][0]*element_coord[0, 2] + diffN[2][1]*element_coord[1, 2] +
        diffN[2][2]*element_coord[2, 2] + diffN[2][3]*element_coord[3, 2] +
        diffN[2][4]*element_coord[4, 2] + diffN[2][5]*element_coord[5, 2] +
        diffN[2][6]*element_coord[6, 2] + diffN[2][7]*element_coord[7, 2]
    )

# ==============================================================================
# FUNÇÕES CDEF AUXILIARES OTIMIZADAS (PARA EVITAR ALOCAÇÕES EM LOOPS)
# ==============================================================================
    
@cdivision(True)
@exceptval(check=False)
@boundscheck(False)
@wraparound(False)         
@nonecheck(False)
cdef inline void C_detJacobi(FLT64 r0, FLT64 r1, FLT64 r2, FLT64 [:, ::1] element_coord, FLT64 *detJ) nogil:
    cdef FLT64 jac[3][3]
    JACOBIANO(r0, r1, r2, element_coord, jac)
    detJ[0] = jac[0][0] * jac[1][1] * jac[2][2] + jac[0][1] * jac[1][2] * jac[2][0] + jac[0][2] * jac[1][0] * jac[2][1] - jac[2][0] * jac[1][1] * jac[0][2] - jac[2][1] * jac[1][2] * jac[0][0] - jac[2][2] * jac[1][0] * jac[0][1]

    
@cdivision(True)
@exceptval(check=False)
@boundscheck(False)
@wraparound(False)         
@nonecheck(False)
cdef inline void C_ShapeFunctions(FLT64 r0, FLT64 r1, FLT64 r2, INT32 nodedof, FLT64[:, ::1] matN) nogil:
    cdef FLT64 N[8]
    MATN(r0, r1, r2, N)
    cdef Py_ssize_t block, dof
    for block in range(8):
        for dof in range(nodedof):
            matN[dof, block * nodedof + dof] = N[block]

    
@cdivision(True)
@exceptval(check=False)
@boundscheck(False)
@wraparound(False)         
@nonecheck(False)
cdef inline void C_DiffShapeFunction(FLT64 r0, FLT64 r1, FLT64 r2, INT32 nodedof, FLT64[:, ::1] matdiffN) nogil:
    cdef FLT64 dN[3][8]
    MATDIFFN(r0, r1, r2, dN)
    cdef Py_ssize_t block, dof
    for block in range(8):
        for dof in range(nodedof):
            # Correção: Uso de colchetes duplos separados para array estático C
            matdiffN[nodedof * dof - dof * (nodedof - 3), block * nodedof + dof] = dN[0][block]
            matdiffN[nodedof * dof - dof * (nodedof - 3) + 1, block * nodedof + dof] = dN[1][block]
            matdiffN[nodedof * dof - dof * (nodedof - 3) + 2, block * nodedof + dof] = dN[2][block]

    
@cdivision(True)
@exceptval(check=False)
@boundscheck(False)
@wraparound(False)         
@nonecheck(False)
cdef inline void C_invJacobi(FLT64 r0, FLT64 r1, FLT64 r2, FLT64[:, ::1] element_coord, INT32 nodedof, FLT64[:, ::1] mat_invJ) nogil:
    cdef FLT64 jac[3][3], invJ[3][3]
    JACOBIANO(r0, r1, r2, element_coord, jac)
    INV(jac, invJ)
    cdef Py_ssize_t block, dimr, dimc
    for block in range(nodedof):
        for dimr in range(3):
            for dimc in range(3):
                # Correção: Acesso correto ao array estático C invJ com [dimr][dimc]
                mat_invJ[block * nodedof + dimr, block * nodedof + dimc] = invJ[dimr][dimc]

    
@cdivision(True)
@exceptval(check=False)
@boundscheck(False)
@wraparound(False)         
@nonecheck(False)
cdef inline void C_compute_B(INT32[:, ::1] H, FLT64[:, ::1] invJ, FLT64[:, ::1] diffN, FLT64[:, ::1] T, FLT64[:, ::1] B) nogil:
    cdef INT32 h_rows = H.shape[0]
    cdef INT32 h_cols = H.shape[1]
    cdef INT32 invJ_cols = invJ.shape[1]
    cdef INT32 diffN_cols = diffN.shape[1]
    cdef Py_ssize_t i, j, k

    for i in range(h_rows):
        for j in range(invJ_cols):
            T[i, j] = 0.0
        for j in range(diffN_cols):
            B[i, j] = 0.0

    for i in range(h_rows):
        for j in range(invJ_cols):
            for k in range(h_cols):
                T[i, j] += H[i, k] * invJ[k, j]

    for i in range(h_rows):
        for j in range(diffN_cols):
            for k in range(invJ_cols):
                B[i, j] += T[i, k] * diffN[k, j]

# ==============================================================================
# FUNÇÕES EXPORTADAS (PYTHON / DEF)
# ==============================================================================
    
@cdivision(True)
@exceptval(check=False)
@boundscheck(False)
@wraparound(False)         
@nonecheck(False)
def ShapeFunctions(FLT64 [::1] point_gauss, INT32 nodedof):
    cdef np.ndarray[FLT64, ndim=2, mode="c"] matN = np.zeros((nodedof, 8 * nodedof), dtype=np.float64) 
    C_ShapeFunctions(point_gauss[0], point_gauss[1], point_gauss[2], nodedof, matN)
    return matN

    
@cdivision(True)
@exceptval(check=False)
@boundscheck(False)
@wraparound(False)         
@nonecheck(False)
def DiffShapeFuntion(FLT64 [::1] point_gauss, INT32 nodedof):
    cdef np.ndarray[FLT64, ndim=2, mode="c"] matdiffN = np.zeros((3 * nodedof, 8 * nodedof), dtype=np.float64) 
    C_DiffShapeFunction(point_gauss[0], point_gauss[1], point_gauss[2], nodedof, matdiffN)
    return matdiffN
    
    
@cdivision(True)
@exceptval(check=False)
@boundscheck(False)
@wraparound(False)         
@nonecheck(False)
def Jacobian(FLT64 [::1] point_gauss, FLT64 [:, ::1] element_coord):  
    cdef FLT64 jac[3][3]
    JACOBIANO(point_gauss[0], point_gauss[1], point_gauss[2], element_coord, jac)
    cdef np.ndarray[FLT64, ndim=2, mode="c"] mat_jac = np.zeros((3, 3), dtype=np.float64)
    mat_jac[0, 0] = jac[0][0]
    mat_jac[0, 1] = jac[0][1]
    mat_jac[0, 2] = jac[0][2]
    mat_jac[1, 0] = jac[1][0]
    mat_jac[1, 1] = jac[1][1]
    mat_jac[1, 2] = jac[1][2]
    mat_jac[2, 0] = jac[2][0]
    mat_jac[2, 1] = jac[2][1]
    mat_jac[2, 2] = jac[2][2]
    return mat_jac
    
    
@cdivision(True)
@exceptval(check=False)
@boundscheck(False)
@wraparound(False)         
@nonecheck(False)
def invJacobi(FLT64 [::1] point_gauss, FLT64 [:, ::1] element_coord, INT32 nodedof):
    cdef np.ndarray[FLT64, ndim=2, mode="c"] mat_invJ = np.zeros((3 * nodedof, 3 * nodedof), dtype=np.float64)  
    C_invJacobi(point_gauss[0], point_gauss[1], point_gauss[2], element_coord, nodedof, mat_invJ)
    return mat_invJ
    
@cdivision(True)
@exceptval(check=False)
@boundscheck(False)
@wraparound(False)         
@nonecheck(False)
def detJacobi(FLT64 [::1] point_gauss, FLT64 [:, ::1] element_coord):
    cdef FLT64 detJ
    C_detJacobi(point_gauss[0], point_gauss[1], point_gauss[2], element_coord, &detJ)
    return detJ

    
@cdivision(True)
@exceptval(check=False)
@boundscheck(False)
@wraparound(False)         
@nonecheck(False)
def compute_VOL(    
    FLT64 [::1] pt_gauss_nodes,
    FLT64 [::1] weight_factor,
    FLT64 [:, ::1] element_coord,
):
    cdef Py_ssize_t ip, jp, kp
    cdef FLT64 detJ, VOL = 0.0
    cdef FLT64 r0, r1, r2

    for ip in range(1):
        for jp in range(1):
            for kp in range(1):
                r0 = pt_gauss_nodes[ip]
                r1 = pt_gauss_nodes[jp]
                r2 = pt_gauss_nodes[kp]
                C_detJacobi(r0, r1, r2, element_coord, &detJ)
                VOL = VOL + fabs(detJ) * weight_factor[ip] * weight_factor[jp] * weight_factor[kp]
    return VOL
    
@cdivision(True)
@exceptval(check=False)
@boundscheck(False)
@wraparound(False)         
@nonecheck(False)
def compute_B(INT32 [:, ::1] H, FLT64 [:, ::1] invJ, FLT64 [:, ::1] diffN):
    cdef INT32 h_rows = H.shape[0]
    cdef INT32 invJ_cols = invJ.shape[1]
    cdef INT32 diffN_cols = diffN.shape[1]
    
    cdef np.ndarray[FLT64, ndim=2, mode="c"] T = np.zeros((h_rows, invJ_cols), dtype=np.float64)
    cdef np.ndarray[FLT64, ndim=2, mode="c"] B = np.zeros((h_rows, diffN_cols), dtype=np.float64)
    
    C_compute_B(H, invJ, diffN, T, B)
    return B

    
@cdivision(True)
@exceptval(check=False)
@boundscheck(False)
@wraparound(False)         
@nonecheck(False)
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
    cdef Py_ssize_t ip, jp, kp
    cdef INT32 r, c, k
    cdef FLT64 detJ, scale, r0, r1, r2
    cdef INT32 c_rows = C.shape[0]

    cdef np.ndarray[FLT64, ndim=2, mode="c"] K_elem_mat = np.zeros((elemdof, elemdof), dtype=np.float64)
    cdef np.ndarray[FLT64, ndim=2, mode="c"] CB = np.zeros((c_rows, elemdof), dtype=np.float64)
    
    cdef np.ndarray[FLT64, ndim=2, mode="c"] diffN = np.zeros((3 * nodedof, 8 * nodedof), dtype=np.float64)
    cdef np.ndarray[FLT64, ndim=2, mode="c"] invJ = np.zeros((3 * nodedof, 3 * nodedof), dtype=np.float64)
    cdef np.ndarray[FLT64, ndim=2, mode="c"] T = np.zeros((H.shape[0], invJ.shape[1]), dtype=np.float64)
    cdef np.ndarray[FLT64, ndim=2, mode="c"] B = np.zeros((H.shape[0], diffN.shape[1]), dtype=np.float64)

    for ip in range(intgauss):
        for jp in range(intgauss):
            for kp in range(intgauss):
                r0 = pt_gauss_nodes[ip]
                r1 = pt_gauss_nodes[jp]
                r2 = pt_gauss_nodes[kp]
                
                C_detJacobi(r0, r1, r2, element_coord, &detJ)
                C_DiffShapeFunction(r0, r1, r2, nodedof, diffN)
                C_invJacobi(r0, r1, r2, element_coord, nodedof, invJ)
                C_compute_B(H, invJ, diffN, T, B)
                
                scale = fabs(detJ) * weight_factor[ip] * weight_factor[jp] * weight_factor[kp]
                
                for r in range(c_rows):
                    for c in range(elemdof):
                        CB[r, c] = 0.0
                        for k in range(c_rows):
                            CB[r, c] += C[r, k] * B[k, c]
                            
                for r in range(elemdof):
                    for c in range(elemdof):
                        for k in range(c_rows):
                            K_elem_mat[r, c] += B[k, r] * CB[k, c] * scale

    return K_elem_mat

    
@cdivision(True)
@exceptval(check=False)
@boundscheck(False)
@wraparound(False)         
@nonecheck(False)
def MassLinear(
    FLT64 [::1] pt_gauss_nodes,      
    FLT64 [::1] weight_factor,       
    INT32 intgauss,
    FLT64 [:, ::1] element_coord,
    INT32 elemdof,
    INT32 nodedof,
    FLT64 R,                       
):
    cdef Py_ssize_t ip, jp, kp
    cdef INT32 r, c, k
    cdef FLT64 detJ, scale, r0, r1, r2
    
    cdef np.ndarray[FLT64, ndim=2, mode="c"] M_elem_mat = np.zeros((elemdof, elemdof), dtype=np.float64)
    cdef np.ndarray[FLT64, ndim=2, mode="c"] N = np.zeros((nodedof, 8 * nodedof), dtype=np.float64)

    for ip in range(intgauss):
        for jp in range(intgauss):
            for kp in range(intgauss):
                r0 = pt_gauss_nodes[ip]
                r1 = pt_gauss_nodes[jp]
                r2 = pt_gauss_nodes[kp]
                
                C_detJacobi(r0, r1, r2, element_coord, &detJ)
                C_ShapeFunctions(r0, r1, r2, nodedof, N) 
                
                scale = R * fabs(detJ) * weight_factor[ip] * weight_factor[jp] * weight_factor[kp]
                
                for r in range(elemdof):
                    for c in range(elemdof):
                        for k in range(nodedof):
                            M_elem_mat[r, c] += N[k, r] * N[k, c] * scale

    return M_elem_mat

    
@cdivision(True)
@exceptval(check=False)
@boundscheck(False)
@wraparound(False)         
@nonecheck(False)
def NodeList(INT32 [:, ::1] inci, INT32 element_number):
    cdef np.ndarray[INT32, ndim=1, mode="c"] node_list = np.zeros(8, dtype=np.int32)
    node_list[0] = inci[element_number, 4]
    node_list[1] = inci[element_number, 5]
    node_list[2] = inci[element_number, 6]
    node_list[3] = inci[element_number, 7]
    node_list[4] = inci[element_number, 8]
    node_list[5] = inci[element_number, 9]
    node_list[6] = inci[element_number, 10]
    node_list[7] = inci[element_number, 11]
    return node_list
            
    
@cdivision(True)
@exceptval(check=False)
@boundscheck(False)
@wraparound(False)         
@nonecheck(False)
def NodeCoord(FLT64 [:, ::1] coord, INT32 [::1] nodelist):
    cdef np.ndarray[FLT64, ndim=2, mode="c"] element_coord = np.zeros((8, 3), dtype=np.float64)
    element_coord[0, 0] = coord[nodelist[0] - 1, 1]
    element_coord[0, 1] = coord[nodelist[0] - 1, 2]
    element_coord[0, 2] = coord[nodelist[0] - 1, 3]
    element_coord[1, 0] = coord[nodelist[1] - 1, 1]
    element_coord[1, 1] = coord[nodelist[1] - 1, 2]
    element_coord[1, 2] = coord[nodelist[1] - 1, 3]
    element_coord[2, 0] = coord[nodelist[2] - 1, 1]
    element_coord[2, 1] = coord[nodelist[2] - 1, 2]
    element_coord[2, 2] = coord[nodelist[2] - 1, 3]
    element_coord[3, 0] = coord[nodelist[3] - 1, 1]
    element_coord[3, 1] = coord[nodelist[3] - 1, 2]
    element_coord[3, 2] = coord[nodelist[3] - 1, 3]
    element_coord[4, 0] = coord[nodelist[4] - 1, 1]
    element_coord[4, 1] = coord[nodelist[4] - 1, 2]
    element_coord[4, 2] = coord[nodelist[4] - 1, 3]
    element_coord[5, 0] = coord[nodelist[5] - 1, 1]
    element_coord[5, 1] = coord[nodelist[5] - 1, 2]
    element_coord[5, 2] = coord[nodelist[5] - 1, 3]
    element_coord[6, 0] = coord[nodelist[6] - 1, 1]
    element_coord[6, 1] = coord[nodelist[6] - 1, 2]
    element_coord[6, 2] = coord[nodelist[6] - 1, 3]
    element_coord[7, 0] = coord[nodelist[7] - 1, 1]
    element_coord[7, 1] = coord[nodelist[7] - 1, 2]
    element_coord[7, 2] = coord[nodelist[7] - 1, 3]
    return element_coord
    
@cdivision(True)
@exceptval(check=False)
@boundscheck(False)
@wraparound(False)         
@nonecheck(False)
def LocKey(INT32 [::1] node_list, INT32 nodedof):
    cdef Py_ssize_t n_nodes = node_list.shape[0]
    cdef np.ndarray[INT32, ndim=1, mode="c"] shape_key = np.zeros(n_nodes * nodedof, dtype=np.int32)
    cdef Py_ssize_t node, dof
    for node in range(n_nodes):
        for dof in range(nodedof):
            shape_key[nodedof * node + dof] = nodedof * node_list[node] - (nodedof - dof)
    return shape_key