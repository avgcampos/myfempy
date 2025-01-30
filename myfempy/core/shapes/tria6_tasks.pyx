cimport openmp
from cython cimport boundscheck, wraparound, cdivision, exceptval, nonecheck
from cython.parallel import parallel, prange

import numpy as np
cimport numpy as np

ctypedef np.int32_t INT32
ctypedef np.float64_t FLT64

@cdivision(True)
@exceptval(check=False)
@boundscheck(False) # turn off bounds-checking for entire function
@wraparound(False)  # turn off negative index wrapping for entire function           
@nonecheck(False) 
cdef FLT64 [:, ::1] MATN(FLT64 [::1] r_coord):
    cdef FLT64 r0 = r_coord[0]
    cdef FLT64 r1 = r_coord[1]
    N = np.zeros((1,6), dtype=np.float64)
    cdef FLT64 [:, ::1] N_view = N
    N_view[0][0] = (1 - r0 - r1)*(1 - 2*r0 - 2*r1)
    N_view[0][1] = r0*(2*r0 - 1)
    N_view[0][2] = r1*(2*r1 - 1)
    N_view[0][3] = 4*r0*(1 - r0 - r1)
    N_view[0][4] = 4*r0*r1
    N_view[0][5] = 4*r1*(1 - r0 - r1)
    return N

@cdivision(True)
@exceptval(check=False)
@boundscheck(False) # turn off bounds-checking for entire function
@wraparound(False)  # turn off negative index wrapping for entire function           
@nonecheck(False) 
cdef FLT64 [:, ::1] MATDIFFN(FLT64 [::1] r):
    cdef FLT64 r0 = r[0]
    cdef FLT64 r1 = r[1]
    dN = np.zeros((2,6), dtype=np.float64)
    cdef FLT64 [:, ::1] dN_view = dN
    dN_view[0, 0] = 4*r0 + 4*r1 - 3
    dN_view[0, 1] = 4*r0 - 1
    dN_view[0, 2] = 0
    dN_view[0, 3] = -8*r0 - 4*r1 + 4
    dN_view[0, 4] = 4*r1
    dN_view[0, 5] = -4*r1
    dN_view[1, 0] = 4*r0 + 4*r1 - 3
    dN_view[1, 1] = 0
    dN_view[1, 2] = 4*r1 - 1
    dN_view[1, 3] = -4*r0
    dN_view[1, 4] = 4*r0
    dN_view[1, 5] = -8*r1 - 4*r0 + 4
    return dN

@cdivision(True)
@exceptval(check=False)
@boundscheck(False) # turn off bounds-checking for entire function
@wraparound(False)  # turn off negative index wrapping for entire function           
@nonecheck(False) 
cdef FLT64 DET(FLT64 [:, ::1] A):
    cdef FLT64 det = A[0][0]*A[1][1]-A[0][1]*A[1][0]
    return det

@cdivision(True)
@exceptval(check=False)
@boundscheck(False) # turn off bounds-checking for entire function
@wraparound(False)  # turn off negative index wrapping for entire function           
@nonecheck(False) 
cdef FLT64 [:, ::1] INV(FLT64 [:, ::1] A):
    cdef FLT64 detA = DET(A)
    cdef FLT64 [:, ::1] invA = (1/detA)*np.array([[A[1][1],-A[0][1]],
                                                  [-A[1][0], A[0][0]]])
    return invA

@cdivision(True)
@exceptval(check=False)
@boundscheck(False) # turn off bounds-checking for entire function
@wraparound(False)  # turn off negative index wrapping for entire function           
@nonecheck(False) 
cdef FLT64 [:, ::1] JACOBIANO(FLT64 [::1] r_coord, FLT64 [:, ::1] element_coord):  
    cdef FLT64 [:, ::1] diffN = MATDIFFN(r_coord)  
    jac = np.zeros((2, 2), dtype=np.float64)
    cdef FLT64 [:, ::1] jac_view = jac
    jac_view[0][0] = diffN[0][0]*element_coord[0][0]+diffN[0][1]*element_coord[1][0]+diffN[0][2]*element_coord[2][0]+diffN[0][3]*element_coord[3][0]+diffN[0][4]*element_coord[4][0]+diffN[0][5]*element_coord[5][0]
    jac_view[0][1] = diffN[0][0]*element_coord[0][1]+diffN[0][1]*element_coord[1][1]+diffN[0][2]*element_coord[2][1]+diffN[0][3]*element_coord[3][1]+diffN[0][4]*element_coord[4][1]+diffN[0][5]*element_coord[5][1]
    jac_view[1][0] = diffN[1][0]*element_coord[0][0]+diffN[1][1]*element_coord[1][0]+diffN[1][2]*element_coord[2][0]+diffN[1][3]*element_coord[3][0]+diffN[1][4]*element_coord[4][0]+diffN[1][5]*element_coord[5][0]
    jac_view[1][1] = diffN[1][0]*element_coord[0][1]+diffN[1][1]*element_coord[1][1]+diffN[1][2]*element_coord[2][1]+diffN[1][3]*element_coord[3][1]+diffN[1][4]*element_coord[4][1]+diffN[1][5]*element_coord[5][1]
    return jac

@boundscheck(False) # turn off bounds-checking for entire function
@wraparound(False)  # turn off negative index wrapping for entire function                  
def ShapeFunctions(FLT64 [::1] r_coord, INT32 nodedof):
    cdef FLT64 [:, ::1] shape_function = MATN(r_coord)
    matN = np.zeros((nodedof, 6*nodedof), dtype=np.float64) 
    cdef FLT64 [:, ::1] matN_view = matN
    cdef Py_ssize_t block, dof
    for block in range(6):
        for dof in range(nodedof):
            matN_view[dof, block*nodedof+dof] = shape_function[0, block]
    return matN

@boundscheck(False) # turn off bounds-checking for entire function
@wraparound(False)  # turn off negative index wrapping for entire function           
def DiffShapeFuntion(FLT64 [::1] r_coord, INT32 nodedof):
    cdef FLT64 [:, ::1] diff_shape_function = MATDIFFN(r_coord)
    matdiffN = np.zeros((2*nodedof, 6*nodedof), dtype=np.float64) 
    cdef FLT64 [:, ::1] matdiffN_view = matdiffN
    cdef Py_ssize_t block, dof
    for block in range(6):
        for dof in range(nodedof):
            matdiffN_view[nodedof*dof-dof*(nodedof-2), block*nodedof+dof] = diff_shape_function[0, block]
            matdiffN_view[nodedof*dof-dof*(nodedof-2)+1, block*nodedof+dof] = diff_shape_function[1, block]
    return matdiffN
    
@boundscheck(False) # turn off bounds-checking for entire function
@wraparound(False)  # turn off negative index wrapping for entire function           
def Jacobian(FLT64 [::1] r_coord, FLT64 [:, ::1] element_coord):  
    cdef FLT64 [:, ::1] Jac = JACOBIANO(r_coord, element_coord)
    return Jac
    
@boundscheck(False) # turn off bounds-checking for entire function
@wraparound(False)  # turn off negative index wrapping for entire function           
def invJacobi(FLT64 [::1] r_coord, FLT64 [:, ::1] element_coord, INT32 nodedof):
    cdef FLT64 [:, ::1] Jac = JACOBIANO(r_coord, element_coord)
    cdef FLT64 [:, ::1] invJ = INV(Jac)
    mat_invJ = np.zeros((2*nodedof, 2*nodedof), dtype=np.float64)  
    cdef FLT64 [:, ::1] mat_invJ_view = mat_invJ
    cdef Py_ssize_t block, dimr, dimc
    for block in range(nodedof):
        for dimr in range(2):
            for dimc in range(2):
                mat_invJ_view[block*nodedof+dimr, block*nodedof+dimc] = invJ[dimr, dimc]
    
    return mat_invJ

@boundscheck(False) # turn off bounds-checking for entire function
@wraparound(False)  # turn off negative index wrapping for entire function           
def detJacobi(FLT64 [::1] r_coord, FLT64 [:, ::1] element_coord):
    cdef FLT64 [:, ::1] Jac = JACOBIANO(r_coord, element_coord)
    cdef FLT64 detJ = 0.5*DET(Jac)
    return detJ

@boundscheck(False) # turn off bounds-checking for entire function
@wraparound(False)  # turn off negative index wrapping for entire function           
def NodeList(INT32 [:, ::1] inci, INT32 element_number):
    cdef INT32 noi = inci[element_number, 4]
    cdef INT32 noj = inci[element_number, 5]
    cdef INT32 nok = inci[element_number, 6]
    cdef INT32 nol = inci[element_number, 7]
    cdef INT32 nom = inci[element_number, 8]
    cdef INT32 non = inci[element_number, 9]
    cdef INT32 [::1] node_list = np.array([noi, noj, nok, nol, nom, non], dtype=np.int32)                  
    return node_list
            
@boundscheck(False) # turn off bounds-checking for entire function
@wraparound(False)  # turn off negative index wrapping for entire function           
def NodeCoord(FLT64 [:, ::1] coord, INT32 [::1] node_list):
    cdef INT32 noi = node_list[0]
    cdef INT32 noj = node_list[1]
    cdef INT32 nok = node_list[2]
    cdef INT32 nol = node_list[3]
    cdef INT32 nom = node_list[4]
    cdef INT32 non = node_list[5]
    cdef FLT64 xi = coord[noi - 1, 1]
    cdef FLT64 yi = coord[noi - 1, 2]
    cdef FLT64 xj = coord[noj - 1, 1]
    cdef FLT64 yj = coord[noj - 1, 2]
    cdef FLT64 xk = coord[nok - 1, 1]
    cdef FLT64 yk = coord[nok - 1, 2]
    cdef FLT64 xl = coord[nol - 1, 1]
    cdef FLT64 yl = coord[nol - 1, 2]
    cdef FLT64 xm = coord[nom - 1, 1]
    cdef FLT64 ym = coord[nom - 1, 2]
    cdef FLT64 xn = coord[non - 1, 1]
    cdef FLT64 yn = coord[non - 1, 2]
    cdef FLT64 [:, ::1] element_coord = np.array([[xi, yi], [xj, yj], [xk, yk], [xl, yl], [xm, ym], [xn, yn]], dtype=np.float64)
    return element_coord

@boundscheck(False) # turn off bounds-checking for entire function
@wraparound(False)  # turn off negative index wrapping for entire function           
def LocKey(INT32 [::1] node_list, INT32 nodedof):
    shape_key = np.zeros(6*nodedof, dtype=np.int32)
    cdef INT32 [::1] shape_key_view = shape_key
    cdef Py_ssize_t node, dof
    for node in range(len(node_list)):
        for dof in range(nodedof):
            shape_key_view[nodedof*node+dof] = nodedof * node_list[node] - (nodedof-dof)
    return shape_key