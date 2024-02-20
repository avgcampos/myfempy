from __future__ import annotations

import numdifftools as nd
import numpy as np

from myfempy.felib.shapes.shape import Shape


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
        # if poly_order == '1':
        N = np.zeros((1,4))
        N[0, 0] = 0.25*(1-r_coord[0])*(1-r_coord[1])
        N[0, 1] = 0.25*(1+r_coord[0])*(1-r_coord[1])
        N[0, 2] = 0.25*(1+r_coord[0])*(1+r_coord[1])
        N[0, 3] = 0.25*(1-r_coord[0])*(1+r_coord[1])
                    
        # elif poly_order == '2':
        #     N = np.zeros((1,8))
        #     N[0, 0] = 0.25*(1-r_coord[0])*(1-r_coord[1])*(-r_coord[0]-r_coord[1]-1)
        #     N[0, 1] = 0.25*(1+r_coord[0])*(1-r_coord[1])*(r_coord[0]-r_coord[1]-1)
        #     N[0, 2] = 0.25*(1+r_coord[0])*(1+r_coord[1])*(+r_coord[0]+r_coord[1]-1)
        #     N[0, 3] = 0.25*(1-r_coord[0])*(1+r_coord[1])*(-r_coord[0]+r_coord[1]-1)
        #     N[0, 4] = 0.5*(1-r_coord[0]**2)*(1-r_coord[1])
        #     N[0, 5] = 0.5*(1+r_coord[0])*(1-r_coord[1]**2)
        #     N[0, 6] = 0.5*(1-r_coord[0]**2)*(1+r_coord[1])
        #     N[0, 7] = 0.5*(1-r_coord[0])*(1-r_coord[1]**2)
        
        return N
        
    def diffN(shape_function, r_coord):
        dN = nd.Gradient(shape_function, n = 1)
        return dN(r_coord)
    
                        
    def getShapeFunctions(r_coord, nodedof):

        shape_function =  Quad4.N(r_coord)
    
        mat_N = np.zeros((nodedof, 4*nodedof)) # np.zeros((2, 8))
        for block in range(4):
            for dof in range(nodedof):
                mat_N[dof, block*nodedof+dof] = shape_function[0, block]
                
        # for dof in range(nodedof):
        #     mat_N[dof, 0*nodedof+dof] = shape_function[0, 0]
        #     mat_N[dof, 1*nodedof+dof] = shape_function[0, 1]
        #     mat_N[dof, 2*nodedof+dof] = shape_function[0, 2]
        #     mat_N[dof, 3*nodedof+dof] = shape_function[0, 3]
                       
        # if nodedof == 1:
        #     mat_N[0, 0] = shape_function[0, 0]
        #     mat_N[0, 1] = shape_function[0, 1]
        #     mat_N[0, 2] = shape_function[0, 2]
        #     mat_N[0, 3] = shape_function[0, 3]
        
        # elif nodedof == 2:
            # mat_N[0, 0] = shape_function[0, 0]
            # mat_N[0, 2] = shape_function[0, 1]
            # mat_N[0, 4] = shape_function[0, 2]
            # mat_N[0, 6] = shape_function[0, 3]
            
            # mat_N[1, 1] = shape_function[0, 0]
            # mat_N[1, 3] = shape_function[0, 1]
            # mat_N[1, 5] = shape_function[0, 2]
            # mat_N[1, 7] = shape_function[0, 3]
        
        # elif nodedof == 3:
        # mat_N = np.zeros((3, 12))
        # mat_N[0, 0] = shape_function[0, 0]
        # mat_N[0, 3] = shape_function[0, 1]
        # mat_N[0, 6] = shape_function[0, 2]
        # mat_N[0, 9] = shape_function[0, 3]
        
        # mat_N[1, 1] = shape_function[0, 0]
        # mat_N[1, 4] = shape_function[0, 1]
        # mat_N[1, 7] = shape_function[0, 2]
        # mat_N[1, 10] = shape_function[0, 3]
        
        # mat_N[2, 2] = shape_function[0, 0]
        # mat_N[2, 5] = shape_function[0, 1]
        # mat_N[2, 8] = shape_function[0, 2]
        # mat_N[2, 11] = shape_function[0, 3]
        
        # else:
        #     pass
        
        return  mat_N

    def getDiffShapeFuntion(shape_function, r_coord, nodedof):
                
        # nodedof = 2

        diff_shape_function = Quad4.diffN(shape_function, r_coord)
        
        mat_diff_N = np.zeros((2*nodedof, 4*nodedof))
        
        for block in range(4):
            for dof in range(nodedof):
                mat_diff_N[nodedof*dof-dof*(nodedof-2), block*nodedof+dof] = diff_shape_function[0, block]
                mat_diff_N[nodedof*dof-dof*(nodedof-2)+1, block*nodedof+dof] = diff_shape_function[1, block]
                
                # mat_d        # for block in range(4):
        #     for dof in range(nodedof):
        #         mat_diff_N[nodedof*dof-dof*(nodedof-2), block*nodedof+dof] = diff_shape_function[0, block]
        #         mat_diff_N[nodedof*dof-dof*(nodedof-2)+1, block*nodedof+dof] = diff_shape_function[1, block]iff_N[dof+1, block*nodedof+dof] = diff_shape_function[1, block]

        # mat_diff_N = np.zeros((4, 12))
             
        # mat_diff_N[0, 1] = diff_shape_function[0, 0]
        # mat_diff_N[1, 2] = diff_shape_function[1, 0]
        # mat_diff_N[2, 1] = diff_shape_function[1, 0]
        # mat_diff_N[2, 2] = diff_shape_function[0, 0]

        # mat_diff_N[0, 4] = diff_shape_function[0, 1]
        # mat_diff_N[1, 5] = diff_shape_function[1, 1]
        # mat_diff_N[2, 4] = diff_shape_function[1, 1]
        # mat_diff_N[2, 5] = diff_shape_function[0, 1]

        # mat_diff_N[0, 7] = diff_shape_function[0, 2]
        # mat_diff_N[1, 8] = diff_shape_function[1, 2]
        # mat_diff_N[2, 7] = diff_shape_function[1, 2]
        # mat_diff_N[2, 8] = diff_shape_function[0, 2]

        # mat_diff_N[0, 10] = diff_shape_function[0, 3]
        # mat_diff_N[1, 11] = diff_shape_function[1, 3]
        # mat_diff_N[2, 10] = diff_shape_function[1, 3]
        # mat_diff_N[2, 11] = diff_shape_function[0, 3]

        # print(mat_diff_N)
        
        return mat_diff_N
    
        
    def getJacobian(shape_function, r_coord, element_coord):    
        # try:
            return np.dot(Quad4.diffN(shape_function, r_coord), element_coord)
        
        # except:
        #     print('Error Line Jacobian Function')   
        
        
    def invJacobi(shape_function, r_coord, element_coord, nodedof):
        
        # nodedof = 2

        J = Quad4.getJacobian(shape_function, r_coord, element_coord)
        
        invJ = np.linalg.inv(J)  
        # mat_invJ = np.kron(np.eye(2, dtype=int), invJ)
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
        J = Quad4.getJacobian(shape_function, r_coord, element_coord)
        return np.linalg.det(J)
    
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
        
        element_coord = np.array([[xi, yi],
                                  [xj, yj],
                                  [xk, yk],
                                  [xl, yl]])
        return element_coord
    
    
    def getShapeKey(node_list, nodedof):
        """element lockey(dof)"""
        
        shape_key = np.zeros(4*nodedof, dtype=int)
        
        for node in range(len(node_list)):
            for dof in range(nodedof):
                shape_key[nodedof*node+dof] = nodedof * node_list[node] - (nodedof-dof)
                
        return shape_key