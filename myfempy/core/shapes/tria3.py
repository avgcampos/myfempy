from __future__ import annotations

import numdifftools as nd
import numpy as np

from myfempy.core.shapes.shape import Shape


class Tria3(Shape):
    '''Triangular 3-Node Shape Class <ConcreteClassService>'''
    
    def getShapeSet():
        shapeset = {
            "def": "3-nodes_conec 1-first_order",
            "key": "tria3",
            "id": 31,
            "nodes": ["i", "j", "k"],
        }
        return shapeset
    
    def N(r_coord):
        N = np.zeros((1,3))
        N[0, 0] = 1 - r_coord[0] - r_coord[1]
        N[0, 1] = r_coord[0]
        N[0, 2] = r_coord[1]
        return N
        
    def diffN(shape_function, r_coord):
        dN = nd.Gradient(shape_function, n = 1)
        return dN(r_coord)
    
                        
    def getShapeFunctions(r_coord, nodedof):
        
        shape_function =  Tria3.N(r_coord)
    
        mat_N = np.zeros((nodedof, 3*nodedof))
        for block in range(3):
            for dof in range(nodedof):
                mat_N[dof, block*nodedof+dof] = shape_function[0, block]
                
        return  mat_N

    def getDiffShapeFuntion(shape_function, r_coord, nodedof):
        
        diff_shape_function = Tria3.diffN(shape_function, r_coord)
        
        mat_diff_N = np.zeros((2*nodedof, 3*nodedof))
        
        for block in range(3):
            for dof in range(nodedof):
                mat_diff_N[nodedof*dof-dof*(nodedof-2), block*nodedof+dof] = diff_shape_function[0, block]
                mat_diff_N[nodedof*dof-dof*(nodedof-2)+1, block*nodedof+dof] = diff_shape_function[1, block]
         
        return mat_diff_N
    
        
    def getJacobian(shape_function, r_coord, element_coord):    
        # try:
            return np.dot(Tria3.diffN(shape_function, r_coord), element_coord)
        
        # except:
        #     print('Error Line Jacobian Function')   
        
        
    def invJacobi(shape_function, r_coord, element_coord, nodedof):
        
        J = Tria3.getJacobian(shape_function, r_coord, element_coord)
        
        invJ = np.linalg.inv(J)
        
        mat_invJ = np.kron(np.eye(nodedof, dtype=int), invJ)
        
        # mat_invJ = np.zeros((4, 4))
        
        # mat_invJ[0, 0] = invJ[0, 0]
        # mat_invJ[0, 1] = invJ[0, 1]
        # mat_invJ[1, 0] = invJ[1, 0]
        # mat_invJ[1, 1] = invJ[1, 1]
        
        # mat_invJ[2, 2] = invJ[0, 0]
        # mat_invJ[2, 3] = invJ[0, 1]
        # mat_invJ[3, 2] = invJ[1, 0]
        # mat_invJ[3, 3] = invJ[1, 1]
        
        return mat_invJ
    
    def detJacobi(shape_function, r_coord, element_coord):
        J = Tria3.getJacobian(shape_function, r_coord, element_coord)
        return np.abs(0.5*np.linalg.det(J))
    
    def getNodeList(inci, element_number):
        
        noi = int(inci[element_number, 4])
        noj = int(inci[element_number, 5])
        nok = int(inci[element_number, 6])
        
        node_list = [noi, noj, nok]          
        
        return node_list
        
    def getNodeCoord(coord, node_list):
        
        noi = node_list[0]
        noj = node_list[1]
        nok = node_list[2]
        
        xi = coord[noi - 1, 1]
        yi = coord[noi - 1, 2]
        xj = coord[noj - 1, 1]
        yj = coord[noj - 1, 2]
        xk = coord[nok - 1, 1]
        yk = coord[nok - 1, 2]
        
        element_coord = np.array([[xi, yi],
                                  [xj, yj],
                                  [xk, yk]])
        return element_coord
    
    
    def getShapeKey(node_list, nodedof):
        """element lockey(dof)"""
        shape_key = np.zeros(3*nodedof, dtype=int)
        
        for node in range(len(node_list)):
            for dof in range(nodedof):
                shape_key[nodedof*node+dof] = nodedof * node_list[node] - (nodedof-dof)
        return shape_key
