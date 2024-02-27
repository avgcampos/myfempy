# from numpy import array, zeros, sqrt, dot, abs, concatenate, int32, float64

import numpy as np
import cython
INT32 = np.uint32
FLT64 = np.float64

from myfempy.core.utilities import gauss_points
from myfempy.core.elements.element import Element

from myfempy.expe.ke_fast.quad4_v2 import Quad4
# from myfempy.core.material.planestress import PlaneStressIsotropic
from myfempy.expe.ke_fast.planestress import PlaneStressIsotropic

class PlaneStressIsoQuad4(Element):
    '''Plane Structural Element Class <ConcreteClassService>'''
                
    def getElementSet():
        
        elemset = {
            "def": "2D-space 2-node_dofs",
            "key": "plane",
            "id": 22,
            "dofs": {'d': {
                        'ux':1,
                        'uy':2},
                     'f': {
                        'fx':1,
                        'fy':2},
            },
            "tensor": ["sxx", "syy", "sxy"],
        }
        return elemset
    

    def getL():
        L = np.array([[1, 0, 0, 0],
                    [0, 0, 0, 1],
                    [0, 1, 1, 0]], dtype=INT32)
        return L

    def getB(Model, elementcoord, ptg, nodedof):
        diffN = Quad4.getDiffShapeFuntion(ptg, nodedof)
        invJ = Quad4.invJacobi(ptg, elementcoord, nodedof)
        invJdiffN = np.dot(invJ, diffN)
        L = PlaneStressIsoQuad4.getL()
        return np.dot(L, invJdiffN)


    def getStifLinearMat(Model, inci, coord, tabmat, tabgeo, intgauss, element_number):

        elem_set = PlaneStressIsoQuad4.getElementSet()
        nodedof = len(elem_set["dofs"]['d'])
        
        shape_set = Quad4.getShapeSet()
        nodecon = len(shape_set['nodes'])
        type_shape = shape_set["key"]    
                
        edof = nodecon * nodedof
        
        nodelist = Quad4.getNodeList(inci, element_number)    
        
        elementcoord = Quad4.getNodeCoord(coord, nodelist)

        E = tabmat[int(inci[element_number, 2]) - 1, 0]  # material elasticity
        v = tabmat[int(inci[element_number, 2]) - 1, 1]  # material poisson ratio
        C = PlaneStressIsotropic.getElasticTensor(E, v)
            
        L = tabgeo[int(inci[element_number, 3] - 1), 4]
            
        pt, wt = gauss_points(type_shape, intgauss)
               
        K_elem_mat = np.zeros((edof, edof), dtype=FLT64)
        for pp in range(intgauss):
            detJ = Quad4.detJacobi(pt[pp], elementcoord)               
            B = PlaneStressIsoQuad4.getB(Model, elementcoord, pt[pp], nodedof) #np.dot(H, np.dot(invJ, diffN))
            BT = B.transpose() #transpose(B)
            BTC = np.dot(BT, C)
            K_elem_mat += np.dot(BTC, B)*L*abs(detJ)*wt[pp]*wt[pp] #multi_dot([BT, C, B])*L*detJ*wt[pp]*wt[pp] #dot(BTC, B)*L*detJ*wt[pp]*wt[pp]     
            # K_elem_mat += (np.transpose(B) @ C @ B)*L*detJ*wt[pp]*wt[pp]
        return K_elem_mat
    
 
    def getMassConsistentMat(Model, inci, coord, tabmat, tabgeo, intgauss, element_number):
        elem_set = PlaneStressIsoQuad4.getElementSet()
        nodedof = len(elem_set["dofs"]['d'])
        shape_set = Quad4.getShapeSet()
        nodecon = len(shape_set['nodes'])
        type_shape = shape_set["key"]    
        edof = nodecon * nodedof
        nodelist = Quad4.getNodeList(inci, element_number)    
        elementcoord = Quad4.getNodeCoord(coord, nodelist)
        R = tabmat[int(inci[element_number, 2]) - 1, 6]  # material density
        L = tabgeo[int(inci[element_number, 3] - 1), 4]
        pt, wt = gauss_points(type_shape, intgauss)
        M_elem_mat = np.zeros((edof, edof),dtype=FLT64)
        for pp in range(intgauss):
            detJ = Quad4.detJacobi(pt[pp], elementcoord)
            N = Quad4.getShapeFunctions(pt[pp], nodedof)
            NT = N.transpose() #transpose(B)
            NTR = np.dot(NT, R)
            M_elem_mat += np.dot(NTR, N)*L*detJ*wt[pp]*wt[pp]  
        return M_elem_mat
    

    def getElementDeformation(U, modelinfo):
        nodetot = modelinfo['nnode']
        nodedof = modelinfo['nodedof']
        Udef = np.zeros((nodetot, 3), dtype=FLT64)
        Umag = np.zeros((nodetot, 1), dtype=FLT64)
        for nn in range(1, nodetot + 1):
            Udef[nn - 1, 0] = U[nodedof * nn - 2]
            Udef[nn - 1, 1] = U[nodedof * nn - 1]
            Umag[nn - 1, 0] = np.sqrt(U[nodedof * nn - 2] ** 2 + U[nodedof * nn - 1] ** 2)
        return np.concatenate((Umag, Udef), axis=1)
        
      
    def setTitleDeformation():
        return ["DISPL_X", "DISPL_Y", "DISPL_Z"] 
    

    def getElementVolume(Model, inci, coord, tabgeo, intgauss, element_number):
        L = tabgeo[int(inci[element_number, 3] - 1), 4]
        shape_set = Quad4.getShapeSet()
        type_shape = shape_set["key"]
        nodelist = Quad4.getNodeList(inci, element_number)
        elementcoord = Quad4.getNodeCoord(coord, nodelist)
        pt, wt = gauss_points(type_shape, intgauss)
        detJ = 0.0
        for pp in range(intgauss):
            N = Quad4.N
            detJ += Quad4.detJacobi(pt[pp], elementcoord)
        return detJ*L