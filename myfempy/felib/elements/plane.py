from __future__ import annotations

import numdifftools as nd
import numpy as np

# from myfempy.felib.material.setpropmat import getElasticity
from myfempy.core.utils import gauss_points
from myfempy.felib.elements.element import Element

# from myfempy.felib.quadrature import gaussian, no_interpol
# from myfempy.felib.materset import get_elasticity


class Plane(Element):
    '''Plane Structural Element Class <ConcreteClassService>'''
    
    # def __init__(self, Shape, Mesh, Material, Geometry) -> None:
    #     self.shape = Shape
    #     self.mesh = Mesh
    #     self.material = Material
    #     self.geometry = Geometry              

    # self.inci =  self.mesh.inci #modelinfo["inci"]
    # self.coord = self.mesh.coord #modelinfo["coord"]
    # self.tabmat = self.material.tabmat #modelinfo["tabmat"]
    # self.tabgeo = self.geometry.tabgeo #modelinfo["tabgeo"]
    # ntensor = self. modelinfo["ntensor"][0]

    # elem_set = Plane.ElementSet()
    # nodedof = len(elem_set["dofs"])

    # shape_set = self.shape.getShapeSet()
    # nodecon = len(shape_set['nodes'])

    # # nodecon = self.shape.nodecon #modelinfo["nodecon"][0]
    # # nodedof = self.shape.nodedof #modelinfo["nodedof"][0]
    # numgaus = self.shape.getNumGaussInt() #modelinfo['intgauss']['ngp']
            
    # # nelem = len(inci)
    # # nnode = len(coord)
    
    # # fulldof = nodedof * nnode
    # edof = nodecon * nodedof
            
    def getElementSet():
        
        elemset = {
            "def": "2D-space 2-node_dofs",
            "key": "plane",
            "id": 22,
            # "dofs": ["ux", "uy"],
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
    
    # def setElementShape(self, Shape):
    #     self.shape = Shape
    
    # def setElementMesh(self, Mesh):
    #     self.mesh = Mesh

    # def setElementMaterial(self, Material):
    #     self.material = Material
        
    # def setElementGeometry(self, Geometry):
    #     self.geometry = Geometry

    def getL():
        L = np.array([[1, 0, 0, 0],
                      [0, 0, 0, 1],
                      [0, 1, 1, 0]])
        return L

    def getB(Model, elementcoord, ptg, nodedof):
        
        diffN = Model.shape.getDiffShapeFuntion(Model.shape.N, ptg, nodedof)
        invJ = Model.shape.invJacobi(Model.shape.N, ptg, elementcoord, nodedof)

        return np.dot(Plane.getL(), np.dot(invJ, diffN))

    def getStifLinearMat(Model, inci, coord, tabmat, tabgeo, intgauss, element_number):

        elem_set = Plane.getElementSet()
        nodedof = len(elem_set["dofs"]['d'])
        
        shape_set = Model.shape.getShapeSet()
        nodecon = len(shape_set['nodes'])
        type_shape = shape_set["key"]    
                
        # nelem = len(inci)
        # nnode = len(coord)
        
        # fulldof = nodedof * nnode
        edof = nodecon * nodedof
        
        nodelist = Model.shape.getNodeList(inci, element_number)    
        
        elementcoord = Model.shape.getNodeCoord(coord, nodelist)

        E = tabmat[int(inci[element_number, 2]) - 1, 0]  # material elasticity
        v = tabmat[int(inci[element_number, 2]) - 1, 1]  # material poisson ratio
        C = Model.material.getElasticTensor(E, v)
            
        L = tabgeo[int(inci[element_number, 3] - 1), 4]
            
        # numgaus = 4 #self.shape.getNumGaussInt()
        pt, wt = gauss_points(type_shape, intgauss)
               
        K_elem_mat = np.zeros((edof, edof))
        for pp in range(intgauss):

            #diffN = Model.shape.getDiffShapeFuntion(Model.shape.N, pt[pp], 1)
            
            detJ = Model.shape.detJacobi(Model.shape.N, pt[pp], elementcoord)
            
            #invJ = Model.shape.invJacobi(Model.shape.N, pt[pp], 1, elementcoord)
                                    
            B = Plane.getB(Model, elementcoord, pt[pp], nodedof) #np.dot(H, np.dot(invJ, diffN))

            K_elem_mat += np.dot(np.dot(np.transpose(B), C), B)*L*detJ*wt[pp]*wt[pp]
        
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
               
        M_elem_mat = np.zeros((edof, edof))
        for pp in range(intgauss):
            detJ = Model.shape.detJacobi(Model.shape.N, pt[pp], elementcoord)
            matN = Model.shape.getShapeFunctions(pt[pp], nodedof)
            M_elem_mat += np.dot(np.dot(np.transpose(matN), R), matN)*L*detJ*wt[pp]*wt[pp]
                     
        return M_elem_mat
    
    def getElementDeformation(U, modelinfo):

        nodetot = modelinfo['nnode']
        nodedof = modelinfo['nodedof']

        Udef = np.zeros((nodetot, 3), dtype=float)
        Umag = np.zeros((nodetot, 1), dtype=float)
        # loading_bar_v1(0, "POST-PROCESSING")
        for nn in range(1, nodetot + 1):
            # loading_bar_v1(100 * ((nn) / self.nnode), "POST-PROCESSING")
            Udef[nn - 1, 0] = U[nodedof * nn - 2]
            Udef[nn - 1, 1] = U[nodedof * nn - 1]
            Umag[nn - 1, 0] = np.sqrt(U[nodedof * nn - 2] ** 2 + U[nodedof * nn - 1] ** 2)
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