# from numpy import array, zeros, eye, dot, matmul, abs, ix_, uint32, float64
import numpy as np
import cython
# from scipy.linalg import inv, kron, det
INT32 = np.uint32
FLT64 = np.float64

# from myfempy.core.utilities import inverse_dim2, determinant_dim2 #, kronProd, determinant, dotProd, getZerosArray, getNewArray, getEyeMatrix
# from myfempy.experimental.utilities import inverse, kronProd, determinant, dotProd, getZerosArray, getNewArray, getEyeMatrix #CYTHON
from myfempy.core.shapes.shape import Shape

class Quad4(Shape):
    '''Quadrilateral 4-Node Shape Class <ConcreteClassService>'''
    
    # @cython.exceptval(check=False)
    # @cython.boundscheck(False)
    # @cython.wraparound(False)
    def getShapeSet():
        shapeset = {
            "def": "4-nodes_conec 1-first_order",
            "key": "quad4",
            "id": 41,
            "nodes": ["i", "j", "k", "l"],
        }
        return shapeset
    
    @cython.exceptval(check=False)
    @cython.boundscheck(False)
    @cython.wraparound(False)    
    def N(r_coord: cython.list):
        N: cython.double[:,:] 
        N = np.zeros((1,4), dtype=FLT64) 
        N[0, 0] = 0.25*(1-r_coord[0])*(1-r_coord[1])
        N[0, 1] = 0.25*(1+r_coord[0])*(1-r_coord[1])
        N[0, 2] = 0.25*(1+r_coord[0])*(1+r_coord[1])
        N[0, 3] = 0.25*(1-r_coord[0])*(1+r_coord[1])
        return N
    
    @cython.exceptval(check=False)
    @cython.boundscheck(False)
    @cython.wraparound(False)   
    def diffN(r: cython.list):
        dN: cython.double[:,:] 
        dN = np.zeros((2,4), dtype=FLT64)
        dN[0, 0] = 0.25*(-1+r[1])
        dN[0, 1] = 0.25*(1-r[1])
        dN[0, 2] = 0.25*(1+r[1])
        dN[0, 3] = 0.25*(-1-r[1])
        dN[1, 0] = 0.25*(-1+r[0])
        dN[1, 1] = 0.25*(-1-r[0])
        dN[1, 2] = 0.25*(1+r[0])
        dN[1, 3] = 0.25*(1-r[0])
        return dN
    
    @cython.exceptval(check=False)
    @cython.boundscheck(False)
    @cython.wraparound(False)                        
    def getShapeFunctions(r_coord: cython.list, nodedof: cython.int):
        shape_function: cython.double[:,:] 
        shape_function = Quad4.N(r_coord)
        mat_N: cython.double[:,:] 
        mat_N = np.zeros((nodedof, 4*nodedof), dtype=FLT64)
        block: cython.int = 0
        dof: cython.int = 0   
        for block in range(4):
            for dof in range(nodedof):
                mat_N[dof, block*nodedof+dof] = shape_function[0, block]
        return  mat_N

    @cython.exceptval(check=False)
    @cython.boundscheck(False)
    @cython.wraparound(False)
    def getDiffShapeFuntion(r_coord: cython.list, nodedof: cython.int):
        diff_shape_function: cython.double[:,:] 
        diff_shape_function = Quad4.diffN(r_coord)
        mat_diff_N: cython.double[:,:] 
        mat_diff_N = np.zeros((2*nodedof, 4*nodedof), dtype=FLT64)
        block: cython.int = 0
        dof: cython.int = 0    
        for block in range(4):
            for dof in range(nodedof):
                mat_diff_N[nodedof*dof-dof*(nodedof-2), block*nodedof+dof] = diff_shape_function[0, block]
                mat_diff_N[nodedof*dof-dof*(nodedof-2)+1, block*nodedof+dof] = diff_shape_function[1, block]
        return mat_diff_N
        
    @cython.exceptval(check=False)
    @cython.boundscheck(False)
    @cython.wraparound(False)        
    def getJacobian(r_coord: cython.list, element_coord: cython.double[:,:]):  
        diffN: cython.double[:,:] 
        diffN = Quad4.diffN(r_coord)  
        jac: cython.double[:,:]
        jac = np.dot(diffN, element_coord) 
        return jac
        
    @cython.exceptval(check=False)
    @cython.boundscheck(False)
    @cython.wraparound(False)
    def invJacobi(r_coord: cython.list, element_coord: cython.double[:,:], nodedof: cython.int):
        J: cython.double[:,:] 
        J = Quad4.getJacobian(r_coord, element_coord)
        invJ: cython.double[:,:] 
        invJ = Quad4.__inverse(np.array(J).flatten())
        mat_invJ: cython.fldoubleoat[:,:] 
        mat_invJ = np.zeros((2*nodedof, 2*nodedof), dtype=FLT64) 
        mat_invJ[0, 0] = invJ[0, 0]
        mat_invJ[0, 1] = invJ[0, 1]
        mat_invJ[1, 0] = invJ[1, 0]
        mat_invJ[1, 1] = invJ[1, 1]
        mat_invJ[2, 2] = invJ[0, 0]
        mat_invJ[2, 3] = invJ[0, 1]
        mat_invJ[3, 2] = invJ[1, 0]
        mat_invJ[3, 3] = invJ[1, 1]     
        return mat_invJ
    
    @cython.exceptval(check=False)
    @cython.boundscheck(False)
    @cython.wraparound(False)    
    def detJacobi(r_coord: cython.list, element_coord: cython.double[:,:]):
        J: cython.fldoubleoat[:,:] 
        J = Quad4.getJacobian(r_coord, element_coord)
        detJ: cython.double 
        detJ = Quad4.__determinant(np.array(J).flatten())
        return detJ
    
    @cython.exceptval(check=False)
    @cython.boundscheck(False)
    @cython.wraparound(False)    
    def getNodeList(inci: cython.double[:,:], element_number: cython.int):
        noi: cython.int = int(inci[element_number, 4])
        noj: cython.int = int(inci[element_number, 5])
        nok: cython.int = int(inci[element_number, 6])
        nol: cython.int = int(inci[element_number, 7])
        node_list: cython.list = [noi, noj, nok, nol]                  
        return node_list
        
    @cython.exceptval(check=False)
    @cython.boundscheck(False)
    @cython.wraparound(False)        
    def getNodeCoord(coord: cython.double[:,:], node_list: cython.list):
        noi: cython.int = node_list[0]
        noj: cython.int = node_list[1]
        nok: cython.int = node_list[2]
        nol: cython.int = node_list[3]
        xi: cython.double = coord[noi - 1, 1]
        yi: cython.double = coord[noi - 1, 2]
        xj: cython.double = coord[noj - 1, 1]
        yj: cython.double = coord[noj - 1, 2]
        xk: cython.double = coord[nok - 1, 1]
        yk: cython.double = coord[nok - 1, 2]
        xl: cython.double = coord[nol - 1, 1]
        yl: cython.double = coord[nol - 1, 2]
        element_coord: cython.double[:,:] = np.array([[xi, yi], [xj, yj], [xk, yk], [xl, yl]], dtype=FLT64)
        return element_coord
    
    @cython.exceptval(check=False)
    @cython.boundscheck(False)
    @cython.wraparound(False)    
    def getShapeKey(node_list: cython.list, nodedof: cython.int):
        """element lockey(dof)"""
        shape_key: cython.ulong[:]
        shape_key = np.zeros(4*nodedof, dtype=INT32)
        node: cython.int = 0
        dof: cython.int = 0  
        for node in range(len(node_list)):
            for dof in range(nodedof):
                shape_key[nodedof*node+dof] = nodedof * node_list[node] - (nodedof-dof)
        return shape_key
    
    @cython.exceptval(check=False)
    @cython.boundscheck(False)
    @cython.wraparound(False)
    def __determinant(A: cython.double[:]):
        detA: cython.double = A[0]*A[3]-A[1]*A[2]
        return detA

    @cython.exceptval(check=False)
    @cython.boundscheck(False)
    @cython.wraparound(False)
    def __inverse(A: cython.double[:]):
        invA: cython.double[:,:] = 1/(A[0]*A[3]-A[1]*A[2])*np.array([[A[3], -A[1]], [-A[2], A[0]]])
        return invA