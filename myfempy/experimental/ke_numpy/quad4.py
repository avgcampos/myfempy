from __future__ import annotations

from numpy import uint32, float64
INT32 = uint32
FLT64 = float64

from myfempy.core.utilities import inverse, kronProd, determinant, dotProd, getZerosArray, getNewArray, getEyeMatrix
# from myfempy.experimental.utilities import inverse, kronProd, determinant, dotProd, getZerosArray, getNewArray, getEyeMatrix #CYTHON
from myfempy.core.shapes.shape import Shape


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
    
    def N(r_coord):
        N = getZerosArray(1, 4, FLT64) #zeros((1,4), dtype=FLT64)
        N[0, 0] = 0.25*(1-r_coord[0])*(1-r_coord[1])
        N[0, 1] = 0.25*(1+r_coord[0])*(1-r_coord[1])
        N[0, 2] = 0.25*(1+r_coord[0])*(1+r_coord[1])
        N[0, 3] = 0.25*(1-r_coord[0])*(1+r_coord[1])

        return N
    
   
    def diffN(N, r):
        dN = getZerosArray(2, 4, FLT64) #zeros((2,4), dtype=FLT64)
        dN[0, 0] = 0.25*(-1+r[1])
        dN[0, 1] = 0.25*(1-r[1])
        dN[0, 2] = 0.25*(1+r[1])
        dN[0, 3] = 0.25*(-1-r[1])
        dN[1, 0] = 0.25*(-1+r[0])
        dN[1, 1] = 0.25*(-1-r[0])
        dN[1, 2] = 0.25*(1+r[0])
        dN[1, 3] = 0.25*(1-r[0])
        return dN
    
                        
    def getShapeFunctions(r_coord, nodedof):
        shape_function =  Quad4.N(r_coord)
        mat_N = getZerosArray(nodedof, 4*nodedof, FLT64) #zeros((nodedof, 4*nodedof), dtype=FLT64)
        for block in range(4):
            for dof in range(nodedof):
                mat_N[dof, block*nodedof+dof] = shape_function[0, block]
        return  mat_N

    def getDiffShapeFuntion(shape_function, r_coord, nodedof):
        diff_shape_function = Quad4.diffN(shape_function, r_coord)
        mat_diff_N = getZerosArray(2*nodedof, 4*nodedof, FLT64) #zeros((2*nodedof, 4*nodedof), dtype=FLT64)
        for block in range(4):
            for dof in range(nodedof):
                mat_diff_N[nodedof*dof-dof*(nodedof-2), block*nodedof+dof] = diff_shape_function[0, block]
                mat_diff_N[nodedof*dof-dof*(nodedof-2)+1, block*nodedof+dof] = diff_shape_function[1, block]
        return mat_diff_N
    
    
    def getJacobian(shape_function, r_coord, element_coord):    
        jac = dotProd(Quad4.diffN(shape_function, r_coord), element_coord) #dot(Quad4.diffN(shape_function, r_coord), element_coord) # comando lento!!!
        return jac
        
    # @profile   
    def invJacobi(shape_function, r_coord, element_coord, nodedof):
        J = Quad4.getJacobian(shape_function, r_coord, element_coord)
        
        # invJ = np.linalg.inv(J)                                 # comando lento!!!
        invJ = inverse(J) #inv(J)
        mat_invJ = getZerosArray(4, 4, FLT64) #zeros((4, 4), dtype=FLT64)  
        # mat_invJ = np.kron(np.eye(nodedof, dtype=int), invJ)    # comando lento!!!
        dim = getEyeMatrix(nodedof, type=INT32)
        mat_invJ = kronProd(dim, invJ) #kron(eye(nodedof, dtype=INT32), invJ)           
        return mat_invJ
    
    def detJacobi(shape_function, r_coord, element_coord):
        J = Quad4.getJacobian(shape_function, r_coord, element_coord)
        detJ = determinant(J)
        return detJ
    
    def getNodeList(inci, element_number):
        
        noi = int(inci[element_number, 4])
        noj = int(inci[element_number, 5])
        nok = int(inci[element_number, 6])
        nol = int(inci[element_number, 7])
        
        node_list = [noi, noj, nok, nol]          
        
        return node_list
        
    def getNodeCoord(coord, node_list):
        noi = node_list[0]
        noj = node_list[1]
        nok = node_list[2]
        nol = node_list[3]
        xi = coord[noi - 1, 1]
        yi = coord[noi - 1, 2]
        xj = coord[noj - 1, 1]
        yj = coord[noj - 1, 2]
        xk = coord[nok - 1, 1]
        yk = coord[nok - 1, 2]
        xl = coord[nol - 1, 1]
        yl = coord[nol - 1, 2]
        array_list = [[xi, yi], [xj, yj], [xk, yk], [xl, yl]]
        element_coord = getNewArray(array_list, FLT64)
        return element_coord
    
    
    def getShapeKey(node_list, nodedof):
        """element lockey(dof)"""
        shape_key = getZerosArray(4*nodedof, 1, FLT64)[:,0].astype(INT32) #zeros(4*nodedof, dtype=INT32)
        for node in range(len(node_list)):
            for dof in range(nodedof):
                shape_key[nodedof*node+dof] = nodedof * node_list[node] - (nodedof-dof)
        return shape_key