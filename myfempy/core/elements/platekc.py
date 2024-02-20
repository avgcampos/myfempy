from __future__ import annotations

import numdifftools as nd
import numpy as np

from myfempy.core.utilities import gauss_points
from myfempy.elements.element import Element

class PlateKC(Element):
    '''Plate Kirchhoff Structural Element Class <ConcreteClassService>'''
                
    def getElementSet():
        
        elemset = {
            "def": "2D-space 3-node_dofs",
            "key": "platekc",
            "id": 23,
            "dofs": {'d': {
                        'uz':1,
                        'rx':2,
                        'ry':3},
                     'f': {
                        'fz':1,
                        'tx':2,
                        'ty':3},
            },
            "tensor": ["sxx", "syy", "sxy"],
        }
        return elemset
    
    def getL():
        # L = np.array([[0, 0, -1, 0, 0, 0],
        #             [0, 0, 0, 0, 0, -1],
        #             [0, 0, 0, -1, -1, 0]])
        
        # derivado w rx ry
        L = np.array([[0, 1, 0, 0, 0, 0],
                      [0, 0, 0, 0, 0, 1],
                      [0, 0, 1, 0, 1, 0]])

        # derivado u v w
        # L = np.array([[-1, 0, 0, 0, 0, 0],
        #               [0, 0, 0, 0, 0, -1],
        #               [0, -1, 0, -1, 0, 0]])
        
        # LIU
        # L = np.array([[0, 0, 0, 0, 1, 0],
        #               [0, 0, 0, 1, 0, 0],
        #               [0, 0, 1, 0, 0, 1]])

        return L

    def getB(Model, elementcoord, ptg, nodedof):
        
        diffN = Model.shape.getDiffShapeFuntion(Model.shape.N, ptg, nodedof)
        invJ = Model.shape.invJacobi(Model.shape.N, ptg, elementcoord, nodedof)

        return np.dot(PlateKC.getL(), np.dot(invJ, diffN))

    def getStifLinearMat(Model, inci, coord, tabmat, tabgeo, intgauss, element_number):

        elem_set = PlateKC.getElementSet()
        nodedof = len(elem_set["dofs"]['d'])

        shape_set = Model.shape.getShapeSet()
        nodecon = len(shape_set['nodes'])
        type_shape = shape_set["key"]    

        edof = nodecon * nodedof
        
        nodelist = Model.shape.getNodeList(inci, element_number)    
        
        elementcoord = Model.shape.getNodeCoord(coord, nodelist)

        E = tabmat[int(inci[element_number, 2]) - 1, 0]  # material elasticity
        v = tabmat[int(inci[element_number, 2]) - 1, 1]  # material poisson ratio
        L = tabgeo[int(inci[element_number, 3] - 1), 4]
            
        C = Model.material.getElasticTensor(E, v)
        C = 0.083333333*L*L*L*C
            
        # numgaus = 4 #self.shape.getNumGaussInt()
        pt, wt = gauss_points(type_shape, intgauss)
               
        K_elem_mat = np.zeros((edof, edof))
        for pp in range(intgauss):
            
            detJ = Model.shape.detJacobi(Model.shape.N, pt[pp], elementcoord)
                                                
            B = PlateKC.getB(Model, elementcoord, pt[pp], nodedof)
            K_elem_mat += np.dot(np.dot(np.transpose(B), C), B)*L*detJ*wt[pp]*wt[pp]
        
        return K_elem_mat
    
    def getMassConsistentMat(Model, inci, coord, tabmat, tabgeo, intgauss, element_number):
        elem_set = PlateKC.getElementSet()
        nodedof = len(elem_set["dofs"]['d'])

        shape_set = Model.shape.getShapeSet()
        nodecon = len(shape_set['nodes'])
        type_shape = shape_set["key"]    

        edof = nodecon * nodedof
        
        nodelist = Model.shape.getNodeList(inci, element_number)    
        
        elementcoord = Model.shape.getNodeCoord(coord, nodelist)
        
        R = tabmat[int(inci[element_number, 2]) - 1, 6]  # material density
        L = tabgeo[int(inci[element_number, 3] - 1), 4]
            
        I = np.zeros((3,3))
        I[0,0] = R*L
        I[1,1] = R*0.083333333*L*L*L
        I[2,2] = R*0.083333333*L*L*L
        
        pt, wt = gauss_points(type_shape, intgauss)
               
        M_elem_mat = np.zeros((edof, edof))
        for pp in range(intgauss):
            detJ = Model.shape.detJacobi(Model.shape.N, pt[pp], elementcoord)
            matN = Model.shape.getShapeFunctions(pt[pp], nodedof)
            M_elem_mat += np.dot(np.dot(np.transpose(matN), I), matN)*L*detJ*wt[pp]*wt[pp]
                     
        return M_elem_mat
    
    def getElementDeformation(U, modelinfo):

        nodetot = modelinfo['nnode']
        nodedof = modelinfo['nodedof']

        Udef = np.zeros((nodetot, 3), dtype=float)
        Umag = np.zeros((nodetot, 1), dtype=float)
        # loading_bar_v1(0, "POST-PROCESSING")
        for nn in range(1, nodetot + 1):
            # loading_bar_v1(100 * ((nn) / self.nnode), "POST-PROCESSING")
            Udef[nn - 1, 0] = U[nodedof * nn - 2] #ux
            Udef[nn - 1, 1] = U[nodedof * nn - 1] #uy
            Udef[nn - 1, 2] = U[nodedof * nn - 3] #uw
            Umag[nn - 1, 0] = np.sqrt(U[nodedof * nn - 3] ** 2 + U[nodedof * nn - 2] ** 2 + U[nodedof * nn - 1] ** 2)
        result = np.concatenate((Umag, Udef), axis=1)
        
        return result
        
    def setTitleDeformation():
        title = ["DISPL_X", "DISPL_Y", "DISPL_Z"]
        return title
    
    def getElementVolume(Model, inci, coord, tabgeo, intgauss, element_number):
        L = tabgeo[int(inci[element_number, 3] - 1), 4]
        shape_set = Model.shape.getShapeSet()
        type_shape = shape_set["key"]
        nodelist = Model.shape.getNodeList(inci, element_number)
        elementcoord = Model.shape.getNodeCoord(coord, nodelist)
        pt, wt = gauss_points(type_shape, intgauss)
        detJ = 0.0
        for pp in range(intgauss):
            detJ += Model.shape.detJacobi(Model.shape.N, pt[pp], elementcoord)
        return detJ*L