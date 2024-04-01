cimport cython
cimport openmp

import numpy as np
cimport numpy as np

INT32 = np.int32
FLT64 = np.float64

ctypedef np.int32_t INT32_t
ctypedef np.float64_t FLT64_t
       
   
cdef double [:, :] MATN(double [:] r_coord):
    cdef double [:, :] N = np.zeros((1,3), dtype=FLT64)
    N[0, 0] = 1 - r_coord[0] - r_coord[1]
    N[0, 1] = r_coord[0]
    N[0, 2] = r_coord[1]
    return N
 
cdef double [:, :] MATDIFFN(double [:] r):
    cdef double r1 = r[1]
    cdef double r0 = r[0]
    cdef double [:, :] dN = np.zeros((2,3), dtype=FLT64)
    dN[0, 0] = -1.0
    dN[0, 1] = 1.0
    dN[1, 0] = -1.0
    dN[1, 2] = 1.0
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
    
    matN = np.zeros((nodedof, 3*nodedof), dtype=FLT64) 
    cdef double [:, :] mat_N = matN
    
    cdef Py_ssize_t block, dof
    
    for block in range(3):
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
    
    matdiffN = np.zeros((2*nodedof, 3*nodedof), dtype=FLT64) 
    cdef double [:, :] mat_diff_N = matdiffN
    
    cdef Py_ssize_t block, dof
    
    for block in range(3):
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

    jac[0,0] = diffN[0,0]*element_coord[0,0]+diffN[0,1]*element_coord[1,0]+diffN[0,2]*element_coord[2,0]
    jac[0,1] = diffN[0,0]*element_coord[0,1]+diffN[0,1]*element_coord[1,1]+diffN[0,2]*element_coord[2,1]
    jac[1,0] = diffN[1,0]*element_coord[0,0]+diffN[1,1]*element_coord[1,0]+diffN[1,2]*element_coord[2,0]
    jac[1,1] = diffN[1,0]*element_coord[0,1]+diffN[1,1]*element_coord[1,1]+diffN[1,2]*element_coord[2,1]
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
   
    mat_invJ[0, 0] = invJ[0, 0]
    mat_invJ[0, 1] = invJ[0, 1]
    mat_invJ[1, 0] = invJ[1, 0]
    mat_invJ[1, 1] = invJ[1, 1]
    mat_invJ[2, 2] = invJ[0, 0]
    mat_invJ[2, 3] = invJ[0, 1]
    mat_invJ[3, 2] = invJ[1, 0]
    mat_invJ[3, 3] = invJ[1, 1]     
    
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
    cdef int noi = int(inci[element_number, 4])
    cdef int noj = int(inci[element_number, 5])
    cdef int nok = int(inci[element_number, 6])
    cdef int [:] node_list = np.array([noi, noj, nok])                  
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
    cdef double xi = coord[noi - 1, 1]
    cdef double yi = coord[noi - 1, 2]
    cdef double xj = coord[noj - 1, 1]
    cdef double yj = coord[noj - 1, 2]
    cdef double xk = coord[nok - 1, 1]
    cdef double yk = coord[nok - 1, 2]
    cdef double [:,:] element_coord = np.array([[xi, yi], [xj, yj], [xk, yk]], dtype=FLT64)
    return element_coord

@cython.cdivision(True)
@cython.exceptval(check=False)
@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function           
@cython.nonecheck(False)     
def LocKey(int [:] node_list, int nodedof):
    cdef int [::1] shape_key = np.zeros(3*nodedof, dtype=INT32)
    cdef Py_ssize_t node, dof
    for node in range(len(node_list)):
        for dof in range(nodedof):
            shape_key[nodedof*node+dof] = nodedof * node_list[node] - (nodedof-dof)
    return shape_key