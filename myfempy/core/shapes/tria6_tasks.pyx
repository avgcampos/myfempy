cimport cython
cimport openmp

import numpy as np
cimport numpy as np

INT32 = np.int32
FLT64 = np.float64

ctypedef np.int32_t INT32_t
ctypedef np.float64_t FLT64_t
          
cdef double [:, :] MATN(double [:] r_coord):
    cdef double r0 = r_coord[0]
    cdef double r1 = r_coord[1]
    cdef double [:, :] N = np.zeros((1,6), dtype=FLT64)
    N[0, 0] = (1 - r0 - r1)*(1 - 2*r0 - 2*r1)
    N[0, 1] = r0*(2*r0 - 1)
    N[0, 2] = r1*(2*r1 - 1)
    N[0, 3] = 4*r0*(1 - r0 - r1)
    N[0, 4] = 4*r0*r1
    N[0, 5] = 4*r1*(1 - r0 - r1)
    return N
 
cdef double [:, :] MATDIFFN(double [:] r):
    cdef double r0 = r[0]
    cdef double r1 = r[1]
    cdef double [:, :] dN = np.zeros((2,6), dtype=FLT64)
    dN[0, 0] = 4*r0 + 4*r1 - 3
    dN[0, 1] = 4*r0 - 1
    dN[0, 2] = 0
    dN[0, 3] = -8*r0 - 4*r1 + 4
    dN[0, 4] = 4*r1
    dN[0, 5] = -4*r1
    dN[1, 0] = 4*r0 + 4*r1 - 3
    dN[1, 1] = 0
    dN[1, 2] = 4*r1 - 1
    dN[1, 3] = -4*r0
    dN[1, 4] = 4*r0
    dN[1, 5] = -8*r1 - 4*r0 + 4
    return dN
    
cdef double DET(double [:] A):
    cdef double det = A[0]*A[3]-A[1]*A[2]
    return det
 
cdef double [:,:] INV(double [:] A):
    cdef double [:,:] invA = (1/(A[0]*A[3]-A[1]*A[2]))*np.array([[A[3], -A[1]], [-A[2], A[0]]])
    return invA
    

@cython.cdivision(True)
@cython.exceptval(check=False)
@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function           
@cython.nonecheck(False)              
def ShapeFunctions(double [:] r_coord, int nodedof):
    
    N = MATN(r_coord)
    cdef double [:, :] shape_function = N
    
    matN = np.zeros((nodedof, 6*nodedof), dtype=FLT64) 
    cdef double [:, :] mat_N = matN
    
    cdef Py_ssize_t block, dof
    
    for block in range(6):
        for dof in range(nodedof):
            mat_N[dof, block*nodedof+dof] = shape_function[0, block]
    return matN

@cython.cdivision(True)
@cython.exceptval(check=False)
@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function           
@cython.nonecheck(False)     
def DiffShapeFuntion(double [:] r_coord, int nodedof):
    
    diffN = MATDIFFN(r_coord)
    cdef double [:, :] diff_shape_function = diffN
    
    matdiffN = np.zeros((2*nodedof, 6*nodedof), dtype=FLT64) 
    cdef double [:, :] mat_diff_N = matdiffN
    
    cdef Py_ssize_t block, dof
    
    for block in range(6):
        for dof in range(nodedof):
            mat_diff_N[nodedof*dof-dof*(nodedof-2), block*nodedof+dof] = diff_shape_function[0, block]
            mat_diff_N[nodedof*dof-dof*(nodedof-2)+1, block*nodedof+dof] = diff_shape_function[1, block]
  
    return matdiffN
    
@cython.cdivision(True)
@cython.exceptval(check=False)
@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function           
@cython.nonecheck(False)     
def Jacobian(double [:] r_coord, double [:, :] element_coord):  
    cdef double [:, :] diffN = MATDIFFN(r_coord)  
    cdef double [:, :] jac = np.zeros((2, 2), dtype=FLT64)

    jac[0,0] = diffN[0,0]*element_coord[0,0]+diffN[0,1]*element_coord[1,0]+diffN[0,2]*element_coord[2,0]+diffN[0,3]*element_coord[3,0]+diffN[0,4]*element_coord[4,0]+diffN[0,5]*element_coord[5,0]
    jac[0,1] = diffN[0,0]*element_coord[0,1]+diffN[0,1]*element_coord[1,1]+diffN[0,2]*element_coord[2,1]+diffN[0,3]*element_coord[3,1]+diffN[0,4]*element_coord[4,1]+diffN[0,5]*element_coord[5,1]
    jac[1,0] = diffN[1,0]*element_coord[0,0]+diffN[1,1]*element_coord[1,0]+diffN[1,2]*element_coord[2,0]+diffN[1,3]*element_coord[3,0]+diffN[1,4]*element_coord[4,0]+diffN[1,5]*element_coord[5,0]
    jac[1,1] = diffN[1,0]*element_coord[0,1]+diffN[1,1]*element_coord[1,1]+diffN[1,2]*element_coord[2,1]+diffN[1,3]*element_coord[3,1]+diffN[1,4]*element_coord[4,1]+diffN[1,5]*element_coord[5,1]
    return jac
    
@cython.cdivision(True)
@cython.exceptval(check=False)
@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function           
@cython.nonecheck(False)     
def invJacobi(double [:] r_coord, double [:, :] element_coord, int nodedof):
    cdef double [:, :] J = Jacobian(r_coord, element_coord)
    cdef double [:, :] invJ = INV(np.array(J).flatten())
    cdef double [:, :] mat_invJ = np.zeros((2*nodedof, 2*nodedof), dtype=FLT64)  

    cdef Py_ssize_t block, dimr, dimc

    for block in range(nodedof):
        for dimr in range(2):
            for dimc in range(2):
                mat_invJ[block*nodedof+dimr, block*nodedof+dimc] = invJ[dimr, dimc]    

    return mat_invJ

@cython.cdivision(True)
@cython.exceptval(check=False)
@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function           
@cython.nonecheck(False)     
def detJacobi(double [:] r_coord, double [:, :] element_coord):
    cdef double [:, :] J = Jacobian(r_coord, element_coord)
    cdef double detJ = 0.0
    detJ = DET(np.array(J).flatten())
    return 0.5*detJ

@cython.cdivision(True)
@cython.exceptval(check=False)
@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function           
@cython.nonecheck(False)     
def NodeList(int [:, :] inci, int element_number):
    cdef int noi = inci[element_number, 4]
    cdef int noj = inci[element_number, 5]
    cdef int nok = inci[element_number, 6]
    cdef int nol = inci[element_number, 7]
    cdef int nom = inci[element_number, 8]
    cdef int non = inci[element_number, 9]
    cdef int [:] node_list = np.array([noi, noj, nok, nol, nom, non])                  
    return node_list
            
@cython.cdivision(True)
@cython.exceptval(check=False)
@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function           
@cython.nonecheck(False)     
def NodeCoord(double [:, :] coord, int [:] node_list):
    cdef int noi = node_list[0]
    cdef int noj = node_list[1]
    cdef int nok = node_list[2]
    cdef int nol = node_list[3]
    cdef int nom = node_list[4]
    cdef int non = node_list[5]
    cdef double xi = coord[noi - 1, 1]
    cdef double yi = coord[noi - 1, 2]
    cdef double xj = coord[noj - 1, 1]
    cdef double yj = coord[noj - 1, 2]
    cdef double xk = coord[nok - 1, 1]
    cdef double yk = coord[nok - 1, 2]
    cdef double xl = coord[nol - 1, 1]
    cdef double yl = coord[nol - 1, 2]
    cdef double xm = coord[nom - 1, 1]
    cdef double ym = coord[nom - 1, 2]
    cdef double xn = coord[non - 1, 1]
    cdef double yn = coord[non - 1, 2]
    cdef double [:,:] element_coord = np.array([[xi, yi], [xj, yj], [xk, yk], [xl, yl], [xm, ym], [xn, yn]], dtype=FLT64)
    return element_coord

@cython.cdivision(True)
@cython.exceptval(check=False)
@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function           
@cython.nonecheck(False)     
def LocKey(int [:] node_list, int nodedof):
    cdef int [::1] shape_key = np.zeros(6*nodedof, dtype=INT32)
    cdef Py_ssize_t node, dof
    for node in range(len(node_list)):
        for dof in range(nodedof):
            shape_key[nodedof*node+dof] = nodedof * node_list[node] - (nodedof-dof)
    return shape_key