from __future__ import annotations

from os import environ
environ['OMP_NUM_THREADS'] = '3'
from numpy import array, zeros, transpose, sqrt, dot, abs, matmul, concatenate, int32, float64
from numpy.linalg import multi_dot
INT32 = int32
FLT64 = float64

from myfempy.core.utilities import gauss_points
from myfempy.core.elements.element import Element

class Plane(Element):
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
        return array([[1, 0, 0, 0],
                    [0, 0, 0, 1],
                    [0, 1, 1, 0]], dtype=INT32)

    # @profile
    def getB(Model, elementcoord, ptg, nodedof):
        # N = Model.shape
        diffN = Model.shape.getDiffShapeFuntion(ptg, nodedof)
        invJ = Model.shape.invJacobi(ptg, elementcoord, nodedof)

        invJdiffN = matmul(invJ, diffN)
        L = Plane.getL()
        return matmul(L, invJdiffN)

    # @profile
    def getStifLinearMat(Model, inci, coord, tabmat, tabgeo, intgauss, element_number):

        elem_set = Plane.getElementSet()
        nodedof = len(elem_set["dofs"]['d'])
        
        shape_set = Model.shape.getShapeSet()
        nodecon = len(shape_set['nodes'])
        type_shape = shape_set["key"]    
                
        edof = nodecon * nodedof
        
        nodelist = Model.shape.getNodeList(inci, element_number)    
        
        elementcoord = Model.shape.getNodeCoord(coord, nodelist)

        E = tabmat[int(inci[element_number, 2]) - 1, 0]  # material elasticity
        v = tabmat[int(inci[element_number, 2]) - 1, 1]  # material poisson ratio
        C = Model.material.getElasticTensor(E, v)
            
        L = tabgeo[int(inci[element_number, 3] - 1), 4]
            
        pt, wt = gauss_points(type_shape, intgauss)
               
        K_elem_mat = zeros((edof, edof), dtype=FLT64)
        for pp in range(intgauss):
            # N = Model.shape.N
            detJ = Model.shape.detJacobi(pt[pp], elementcoord)               
            B = Plane.getB(Model, elementcoord, pt[pp], nodedof) #np.dot(H, np.dot(invJ, diffN))
            BT = B.transpose() #transpose(B)
            BTC = dot(BT, C)
            K_elem_mat += dot(BTC, B)*L*abs(detJ)*wt[pp]*wt[pp] #multi_dot([BT, C, B])*L*detJ*wt[pp]*wt[pp] #dot(BTC, B)*L*detJ*wt[pp]*wt[pp]     
            # K_elem_mat += (np.transpose(B) @ C @ B)*L*detJ*wt[pp]*wt[pp]
        return K_elem_mat
    
    def getMassConsistentMat(Model, inci, coord, tabmat, tabgeo, intgauss, element_number):
        elem_set = Plane.getElementSet()
        nodedof = len(elem_set["dofs"]['d'])
        shape_set = Model.shape.getShapeSet()
        nodecon = len(shape_set['nodes'])
        type_shape = shape_set["key"]    
        edof = nodecon * nodedof
        nodelist = Model.shape.getNodeList(inci, element_number)    
        elementcoord = Model.shape.getNodeCoord(coord, nodelist)
        R = tabmat[int(inci[element_number, 2]) - 1, 6]  # material density
        L = tabgeo[int(inci[element_number, 3] - 1), 4]
        pt, wt = gauss_points(type_shape, intgauss)
        M_elem_mat = zeros((edof, edof),dtype=FLT64)
        for pp in range(intgauss):
            # N = Model.shape.N
            detJ = Model.shape.detJacobi(pt[pp], elementcoord)
            matN = Model.shape.getShapeFunctions(pt[pp], nodedof)
            M_elem_mat += multi_dot([transpose(matN), R, matN])*L*abs(detJ)*wt[pp]*wt[pp] #dot(dot(transpose(matN), R), matN)*L*detJ*wt[pp]*wt[pp]  
        return M_elem_mat
    
    def getElementDeformation(U, modelinfo):
        nodetot = modelinfo['nnode']
        nodedof = modelinfo['nodedof']
        Udef = zeros((nodetot, 3), dtype=FLT64)
        Umag = zeros((nodetot, 1), dtype=FLT64)
        for nn in range(1, nodetot + 1):
            Udef[nn - 1, 0] = U[nodedof * nn - 2]
            Udef[nn - 1, 1] = U[nodedof * nn - 1]
            Umag[nn - 1, 0] = sqrt(U[nodedof * nn - 2] ** 2 + U[nodedof * nn - 1] ** 2)
        return concatenate((Umag, Udef), axis=1)
        
    def setTitleDeformation():
        return ["DISPL_X", "DISPL_Y", "DISPL_Z"] 
    
    def getElementVolume(Model, inci, coord, tabgeo, intgauss, element_number):
        L = tabgeo[int(inci[element_number, 3] - 1), 4]
        shape_set = Model.shape.getShapeSet()
        type_shape = shape_set["key"]
        nodelist = Model.shape.getNodeList(inci, element_number)
        elementcoord = Model.shape.getNodeCoord(coord, nodelist)
        pt, wt = gauss_points(type_shape, intgauss)
        detJ = 0.0
        for pp in range(intgauss):
            N = Model.shape.N
            detJ += Model.shape.detJacobi(pt[pp], elementcoord)
        return detJ*L