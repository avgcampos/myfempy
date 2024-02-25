import numpy as np
cimport numpy as np
from cython.parallel cimport prange
cimport cython

INT32 = np.uint32
FLT64 = np.float64

from myfempy.core.shapes.shape import Shape
from myfempy.core.utilities import inverse_dim2, determinant_dim2

class Quad4(Shape):
    '''Quadrilateral 4-Node Shape Class <ConcreteClassService>'''
    
    def getShapeSet():
        shapeset = {
            "def": "4-nodes_conec 1-first_order",
            "key": "quad4",
            "id": 41,
            "nodes": ["i", "j", "k", "l"],
        }
        return shapeset
    
    @cython.exceptval(check=False)
    @cython.boundscheck(False) # turn off bounds-checking for entire function
    @cython.wraparound(False)  # turn off negative index wrapping for entire function
    def N(list r_coord):
        cdef double [:,:] N = np.zeros((1,4), dtype=FLT64)
        N[0, 0] = 0.25*(1.0-r_coord[0])*(1.0-r_coord[1])
        N[0, 1] = 0.25*(1.0+r_coord[0])*(1.0-r_coord[1])
        N[0, 2] = 0.25*(1.0+r_coord[0])*(1.0+r_coord[1])
        N[0, 3] = 0.25*(1.0-r_coord[0])*(1.0+r_coord[1])
        return np.array(N)
    
    
    @cython.exceptval(check=False)
    @cython.boundscheck(False) # turn off bounds-checking for entire function
    @cython.wraparound(False)  # turn off negative index wrapping for entire function   
    def diffN(list r):
        cdef double [:,:] dN = np.zeros((2,4), dtype=FLT64)
        dN[0, 0] = 0.25*(-1.0+r[1])
        dN[0, 1] = 0.25*(1.0-r[1])
        dN[0, 2] = 0.25*(1.0+r[1])
        dN[0, 3] = 0.25*(-1.0-r[1])
        dN[1, 0] = 0.25*(-1.0+r[0])
        dN[1, 1] = 0.25*(-1.0-r[0])
        dN[1, 2] = 0.25*(1.0+r[0])
        dN[1, 3] = 0.25*(1.0-r[0])
        return np.array(dN)
    
    @cython.exceptval(check=False)
    @cython.boundscheck(False) # turn off bounds-checking for entire function
    @cython.wraparound(False)  # turn off negative index wrapping for entire function                        
    def getShapeFunctions(list r_coord, int nodedof):
        cdef double [:,:] shape_function = Quad4.N(r_coord)
        cdef double [:,:] mat_N = np.zeros((nodedof, 4*nodedof), dtype=FLT64) 
        cdef int block, dof
        for block in range(4):
            for dof in range(nodedof):
                mat_N[dof, block*nodedof+dof] = shape_function[0, block]
        return np.array(mat_N)

    @cython.exceptval(check=False)
    @cython.boundscheck(False) # turn off bounds-checking for entire function
    @cython.wraparound(False)  # turn off negative index wrapping for entire function
    def getDiffShapeFuntion(list r_coord, int nodedof):
        cdef double [:,:] diff_shape_function = Quad4.diffN(r_coord)
        cdef double [:,:] mat_diff_N = np.zeros((2*nodedof, 4*nodedof), dtype=FLT64) 
        cdef int block, dof
        for block in range(4):
            for dof in range(nodedof):
                mat_diff_N[nodedof*dof-dof*(nodedof-2), block*nodedof+dof] = diff_shape_function[0, block]
                mat_diff_N[nodedof*dof-dof*(nodedof-2)+1, block*nodedof+dof] = diff_shape_function[1, block]
        return np.array(mat_diff_N)
        
    @cython.exceptval(check=False)
    @cython.boundscheck(False) # turn off bounds-checking for entire function
    @cython.wraparound(False)  # turn off negative index wrapping for entire function        
    def getJacobian(list r_coord, double [:, :] element_coord):  
        cdef double [:,:] diffN = Quad4.diffN(r_coord)  
        cdef double [:,:] jac = np.zeros((2, 2), dtype=FLT64)
        jac = np.dot(diffN, element_coord) #matmul(diffN, element_coord)
        return np.array(jac)
        

    @cython.exceptval(check=False)
    @cython.boundscheck(False) # turn off bounds-checking for entire function
    @cython.wraparound(False)  # turn off negative index wrapping for entire function
    def invJacobi(list r_coord, double [:,:] element_coord, int nodedof):
        cdef double [:,:] J = Quad4.getJacobian(r_coord, element_coord)
        cdef double [:,:] invJ = np.zeros((nodedof, nodedof), dtype=FLT64)
        cdef unsigned long [:,:] dim = np.eye(nodedof, dtype=INT32) 
        cdef double [:,:] mat_invJ = np.zeros((4, 4), dtype=FLT64) 

        invJ = inverse_dim2(np.array(J).flatten())
        
        mat_invJ[0, 0] = invJ[0, 0]
        mat_invJ[0, 1] = invJ[0, 1]
        mat_invJ[1, 0] = invJ[1, 0]
        mat_invJ[1, 1] = invJ[1, 1]
        mat_invJ[2, 2] = invJ[0, 0]
        mat_invJ[2, 3] = invJ[0, 1]
        mat_invJ[3, 2] = invJ[1, 0]
        mat_invJ[3, 3] = invJ[1, 1]         
        return np.array(mat_invJ)
    
    @cython.exceptval(check=False)
    @cython.boundscheck(False) # turn off bounds-checking for entire function
    @cython.wraparound(False)  # turn off negative index wrapping for entire function    
    def detJacobi(list r_coord, double [:,:] element_coord):
        cdef double [:,:] J = Quad4.getJacobian(r_coord, element_coord)
        cdef float detJ = 0.0
        detJ = determinant_dim2(np.array(J).flatten())
        return np.array(detJ)
    
    @cython.exceptval(check=False)
    @cython.boundscheck(False) # turn off bounds-checking for entire function
    @cython.wraparound(False)  # turn off negative index wrapping for entire function    
    def getNodeList(double [:,:] inci, int element_number):
        cdef int noi = int(inci[element_number, 4])
        cdef int noj = int(inci[element_number, 5])
        cdef int nok = int(inci[element_number, 6])
        cdef int nol = int(inci[element_number, 7])
        cdef list node_list = [noi, noj, nok, nol]          
        return node_list
        
    @cython.exceptval(check=False)
    @cython.boundscheck(False) # turn off bounds-checking for entire function
    @cython.wraparound(False)  # turn off negative index wrapping for entire function        
    def getNodeCoord(double [:,:] coord, list node_list):
        cdef int noi = node_list[0]
        cdef int noj = node_list[1]
        cdef int nok = node_list[2]
        cdef int nol = node_list[3]
        cdef float xi = coord[noi - 1, 1]
        cdef float yi = coord[noi - 1, 2]
        cdef float xj = coord[noj - 1, 1]
        cdef float yj = coord[noj - 1, 2]
        cdef float xk = coord[nok - 1, 1]
        cdef float yk = coord[nok - 1, 2]
        cdef float xl = coord[nol - 1, 1]
        cdef float yl = coord[nol - 1, 2]
        cdef double [:,:] element_coord = np.array([[xi, yi], [xj, yj], [xk, yk], [xl, yl]], dtype=FLT64)
        return np.array(element_coord)
    
    @cython.exceptval(check=False)
    @cython.boundscheck(False) # turn off bounds-checking for entire function
    @cython.wraparound(False)  # turn off negative index wrapping for entire function    
    def getShapeKey(list node_list, int nodedof):
        """element lockey(dof)"""
        cdef unsigned long [:] shape_key = np.zeros(4*nodedof, dtype=INT32)
        cdef int node, dof
        for node in range(len(node_list)):
            for dof in range(nodedof):
                shape_key[nodedof*node+dof] = nodedof * node_list[node] - (nodedof-dof)
        return np.array(shape_key)
