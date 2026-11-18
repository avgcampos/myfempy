from cython cimport boundscheck, cdivision, exceptval, nonecheck, wraparound
from libc.math cimport fabs, sqrt

import numpy as np
cimport numpy as np

np.import_array()

ctypedef np.int32_t INT32
ctypedef np.float64_t FLT64

# ==============================================================================
# ELEMENTO DE LINHA 2 NÓS: MOLA 3D (AXIAL + TORÇÃO) + CORDA (TRANSVERSAL)
#
# GDL por nó (nodedof = 4), no sistema LOCAL do elemento (x local = eixo da linha):
#     [ u , v , w , rx ]
#       u  -> axial            : mola linear      (rigidez k_ax)
#       v  -> transversal y    : corda            (tração T)
#       w  -> transversal z    : corda            (tração T)
#       rx -> torção em torno de x : mola de torção (rigidez k_tor)
#
# GDL do elemento (elemdof = 8):
#     [u1 v1 w1 rx1 | u2 v2 w2 rx2]
#
# Campos (linhas de N, B, C, R) -> 0:axial, 1:v, 2:w, 3:torção.
# Todas as interpolações são LINEARES (C0): não há GDL de rotação de flexão
# (a corda não tem rigidez à flexão).
#
# Matrizes constitutivas esperadas (mesmo papel do C e R do frame original):
#     C = diag( k_ax*L , T , T , k_tor*L )     -> K_ax = k_ax , K_corda = T/L , K_tor = k_tor
#     R = diag( mu , mu , mu , iota )          -> mu = massa/comprimento,
#                                                 iota = inércia torcional/comprimento
# ==============================================================================

# ==============================================================================
# FUNÇÕES GEOMÉTRICAS E MATEMÁTICAS CDEF (INTERNAS)
# ==============================================================================

@cdivision(True)
@exceptval(check=False)
@boundscheck(False)
@wraparound(False)
@nonecheck(False)
cdef inline np.ndarray[FLT64, ndim=2, mode="c"] MATN(FLT64 r0):
    cdef np.ndarray[FLT64, ndim=2, mode="c"] N = np.zeros((4, 8), dtype=np.float64)
    cdef FLT64 n1 = 0.5 * (1.0 - r0)
    cdef FLT64 n2 = 0.5 * (1.0 + r0)
    N[0, 0] = n1
    N[0, 4] = n2
    N[1, 1] = n1
    N[1, 5] = n2
    N[2, 2] = n1
    N[2, 6] = n2
    N[3, 3] = n1
    N[3, 7] = n2
    return N

@cdivision(True)
@exceptval(check=False)
@boundscheck(False)
@wraparound(False)
@nonecheck(False)
cdef inline np.ndarray[FLT64, ndim=2, mode="c"] MATDIFFN(FLT64 r0):
    # derivada da geometria (x,y,z) -> usada só no Jacobiano (igual ao original)
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
    # matriz "operador" dos campos: d/dr de cada campo (só 1a derivada, em todas
    # as linhas). Substitui a matriz com 2a derivada de flexão do frame original.
    cdef np.ndarray[FLT64, ndim=2, mode="c"] ddN = np.zeros((4, 8), dtype=np.float64)
    ddN[0, 0] = -0.5
    ddN[0, 4] = 0.5
    ddN[1, 1] = -0.5
    ddN[1, 5] = 0.5
    ddN[2, 2] = -0.5
    ddN[2, 6] = 0.5
    ddN[3, 3] = -0.5
    ddN[3, 7] = 0.5
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
    # J = dx/dr = L/2, com L = norma 3D de (x2 - x1).
    # (o original usava só a coordenada x; aqui funciona para qualquer orientação)
    cdef FLT64 dx = element_coord[3, 0] - element_coord[0, 0]
    cdef FLT64 dy = element_coord[4, 0] - element_coord[1, 0]
    cdef FLT64 dz = element_coord[5, 0] - element_coord[2, 0]
    cdef np.ndarray[FLT64, ndim=2, mode="c"] jac = np.zeros((1, 1), dtype=np.float64)
    jac[0, 0] = 0.5 * sqrt(dx * dx + dy * dy + dz * dz)
    return jac

# ==============================================================================
# FUNÇÕES EXPORTADAS (PYTHON / DEF)
# ==============================================================================

@cdivision(True)
@exceptval(check=False)
@boundscheck(False)
@wraparound(False)
@nonecheck(False)
def ShapeFunctions(FLT64 [::1] r_coord, INT32 nodedof, FLT64 detJ):
    # detJ mantido só por compatibilidade de assinatura (não há Hermite para escalar)
    return MATN(r_coord[0])

@cdivision(True)
@exceptval(check=False)
@boundscheck(False)
@wraparound(False)
@nonecheck(False)
def DiffShapeFuntion(FLT64 [::1] r_coord, INT32 nodedof):
    return MATDIFFN(r_coord[0])

@cdivision(True)
@exceptval(check=False)
@boundscheck(False)
@wraparound(False)
@nonecheck(False)
def DiffDiffShapeFuntion(FLT64 [::1] r_coord, INT32 nodedof, FLT64 detJ):
    return MATDIFFDIFFN(r_coord[0])

@cdivision(True)
@exceptval(check=False)
@boundscheck(False)
@wraparound(False)
@nonecheck(False)
def Jacobian(FLT64 [::1] r_coord, FLT64 [:, ::1] element_coord):
    return JACOBIANO(r_coord[0], element_coord)

@cdivision(True)
@exceptval(check=False)
@boundscheck(False)
@wraparound(False)
@nonecheck(False)
def invJacobi(FLT64 [::1] r_coord, FLT64 [:, ::1] element_coord, INT32 nodedof):
    # todos os campos são de 1a derivada -> inv(J) em toda a diagonal
    cdef np.ndarray[FLT64, ndim=2, mode="c"] Jac = JACOBIANO(r_coord[0], element_coord)
    cdef FLT64 invJ = INV(Jac)
    cdef np.ndarray[FLT64, ndim=2, mode="c"] mat_invJ = np.zeros((4, 4), dtype=np.float64)
    mat_invJ[0, 0] = invJ
    mat_invJ[1, 1] = invJ
    mat_invJ[2, 2] = invJ
    mat_invJ[3, 3] = invJ
    return mat_invJ

@cdivision(True)
@exceptval(check=False)
@boundscheck(False)
@wraparound(False)
@nonecheck(False)
def detJacobi(FLT64 [::1] r_coord, FLT64 [:, ::1] element_coord):
    cdef np.ndarray[FLT64, ndim=2, mode="c"] Jac = JACOBIANO(r_coord[0], element_coord)
    return Jac[0, 0]

@cdivision(True)
@exceptval(check=False)
@boundscheck(False)
@wraparound(False)
@nonecheck(False)
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
    # K = integral( B^T C B ) dx   com   B = H * inv(J) * dN/dr
    # Com C = diag(k_ax*L, T, T, k_tor*L)  ->  K = k * [[1,-1],[-1,1]] por campo
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
    FLT64 [:, ::1] R,
):
    # M = integral( N^T R N ) dx
    #  - pontos de Gauss (2 pts)        -> massa CONSISTENTE  mu*L/6*[[2,1],[1,2]]
    #  - pontos de Lobatto (-1,+1; w=1) -> massa CONCENTRADA  mu*L/2*[[1,0],[0,1]]
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

@cdivision(True)
@exceptval(check=False)
@boundscheck(False)
@wraparound(False)
@nonecheck(False)
def RotationMatrix(FLT64 [:, ::1] element_coord):
    # T (8x8) tal que  u_local = T @ u_global  e  K_global = T^T K_local T.
    # Linhas de T3 = [ex; ey; ez]. Como a corda é isotrópica na transversal
    # (C[1,1] = C[2,2], R[1,1] = R[2,2]), a escolha de ey/ez (perpendiculares a ex)
    # é arbitrária. O GDL de torção rx é escalar em torno de ex (cópia direta).
    cdef FLT64 dx = element_coord[3, 0] - element_coord[0, 0]
    cdef FLT64 dy = element_coord[4, 0] - element_coord[1, 0]
    cdef FLT64 dz = element_coord[5, 0] - element_coord[2, 0]
    cdef FLT64 L = sqrt(dx * dx + dy * dy + dz * dz)
    cdef FLT64 ex0 = dx / L
    cdef FLT64 ex1 = dy / L
    cdef FLT64 ex2 = dz / L
    cdef FLT64 ey0, ey1, ey2, ez0, ez1, ez2, nrm
    cdef np.ndarray[FLT64, ndim=2, mode="c"] T = np.zeros((8, 8), dtype=np.float64)
    cdef Py_ssize_t b, off

    # vetor de referência: z global (ou y global se o elemento for quase vertical)
    if fabs(ex2) < 0.9:
        # ey = z_global x ex
        ey0 = -ex1
        ey1 = ex0
        ey2 = 0.0
    else:
        # ey = y_global x ex
        ey0 = ex2
        ey1 = 0.0
        ey2 = -ex0
    nrm = sqrt(ey0 * ey0 + ey1 * ey1 + ey2 * ey2)
    ey0 = ey0 / nrm
    ey1 = ey1 / nrm
    ey2 = ey2 / nrm

    # ez = ex x ey
    ez0 = ex1 * ey2 - ex2 * ey1
    ez1 = ex2 * ey0 - ex0 * ey2
    ez2 = ex0 * ey1 - ex1 * ey0

    for b in range(2):
        off = 4 * b
        T[off, off]         = ex0
        T[off, off + 1]     = ex1
        T[off, off + 2]     = ex2
        T[off + 1, off]     = ey0
        T[off + 1, off + 1] = ey1
        T[off + 1, off + 2] = ey2
        T[off + 2, off]     = ez0
        T[off + 2, off + 1] = ez1
        T[off + 2, off + 2] = ez2
        T[off + 3, off + 3] = 1.0

    return T

@cdivision(True)
@exceptval(check=False)
@boundscheck(False)
@wraparound(False)
@nonecheck(False)
def NodeList(INT32 [:, ::1] inci, INT32 element_number):
    cdef np.ndarray[INT32, ndim=1, mode="c"] node_list = np.zeros(2, dtype=np.int32)
    node_list[0] = inci[element_number, 4]
    node_list[1] = inci[element_number, 5]
    return node_list

@cdivision(True)
@exceptval(check=False)
@boundscheck(False)
@wraparound(False)
@nonecheck(False)
def NodeCoord(FLT64 [:, ::1] coord, INT32 [::1] node_list):
    cdef np.ndarray[FLT64, ndim=2, mode="c"] element_coord = np.zeros((6, 1), dtype=np.float64)
    element_coord[0, 0] = coord[node_list[0] - 1, 1]
    element_coord[1, 0] = coord[node_list[0] - 1, 2]
    element_coord[2, 0] = coord[node_list[0] - 1, 3]
    element_coord[3, 0] = coord[node_list[1] - 1, 1]
    element_coord[4, 0] = coord[node_list[1] - 1, 2]
    element_coord[5, 0] = coord[node_list[1] - 1, 3]
    return element_coord

@cdivision(True)
@exceptval(check=False)
@boundscheck(False)
@wraparound(False)
@nonecheck(False)
def LocKey(INT32 [::1] node_list, INT32 nodedof):
    cdef np.ndarray[INT32, ndim=1, mode="c"] shape_key = np.zeros(2 * nodedof, dtype=np.int32)
    cdef Py_ssize_t node, dof
    for node in range(2):
        for dof in range(nodedof):
            shape_key[nodedof * node + dof] = nodedof * node_list[node] - (nodedof - dof)
    return shape_key
