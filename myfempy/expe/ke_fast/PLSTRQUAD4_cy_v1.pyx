import numpy as np
cimport numpy as np
# from cython.parallel cimport prange
cimport cython

INT32 = np.uint32
FLT64 = np.float64

from myfempy.core.utilities import gauss_points
from myfempy.core.elements.element import Element

from myfempy.expe.ke_fast.quad4_v2 import Quad4
from myfempy.core.material.planestress import PlaneStressIsotropic

class PlaneStressIsoQuad4(Element):
    '''Plane Structural Element Class <ConcreteClassService>'''
                
    @cython.exceptval(check=False)
    @cython.boundscheck(False) # turn off bounds-checking for entire function
    @cython.wraparound(False)  # turn off negative index wrapping for entire function            
    def getElementSet():
        cdef dict elemset = {
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
    
    @cython.exceptval(check=False)
    @cython.boundscheck(False) # turn off bounds-checking for entire function
    @cython.wraparound(False)  # turn off negative index wrapping for entire function 
    def getL():
        cdef unsigned int [:,:] L = np.array([[1, 0, 0, 0],
                                    [0, 0, 0, 1],
                                    [0, 1, 1, 0]], dtype=INT32)
        return L

    @cython.exceptval(check=False)
    @cython.boundscheck(False) # turn off bounds-checking for entire function
    @cython.wraparound(False)  # turn off negative index wrapping for entire function 
    def getB(Model, double [:,:] elementcoord, list ptg, int nodedof):
        cdef double [:,:] diffN = Quad4.getDiffShapeFuntion(ptg, nodedof)
        cdef double [:,:] invJ = Quad4.invJacobi(ptg, elementcoord, nodedof)
        cdef double [:,:] invJdiffN = np.dot(invJ, diffN)
        cdef unsigned int [:,:] L = PlaneStressIsoQuad4.getL()
        cdef double [:,:] B = np.dot(L, invJdiffN)
        return B

    @cython.exceptval(check=False)
    @cython.boundscheck(False) # turn off bounds-checking for entire function
    @cython.wraparound(False)  # turn off negative index wrapping for entire function 
    def getStifLinearMat(Model, double [:,:] inci, double [:,:] coord, double [:,:] tabmat, double [:,:] tabgeo, int intgauss, int element_number):
        cdef dict elem_set = PlaneStressIsoQuad4.getElementSet()
        cdef unsigned int nodedof = len(elem_set["dofs"]['d'])
        cdef dict shape_set = Quad4.getShapeSet()
        cdef unsigned int nodecon = len(shape_set['nodes'])
        cdef str type_shape = shape_set["key"]     
        cdef unsigned int edof = nodecon * nodedof
        cdef list nodelist = Quad4.getNodeList(inci, element_number)    
        cdef double [:,:] elementcoord = Quad4.getNodeCoord(coord, nodelist)
        cdef float E = tabmat[int(inci[element_number, 2]) - 1, 0]  # material elasticity
        cdef float v = tabmat[int(inci[element_number, 2]) - 1, 1]  # material poisson ratio
        cdef double [:,:] C = PlaneStressIsotropic.getElasticTensor(E, v)
        cdef float L = tabgeo[int(inci[element_number, 3] - 1), 4]
        cdef list pt
        cdef list wt
        pt, wt = gauss_points(type_shape, intgauss)
        cdef double [:,:] K_elem_mat = np.zeros((edof, edof), dtype=FLT64)
        cdef float detJ
        cdef double [:,:] B
        cdef double [:,:] BT
        cdef double [:,:] BTC
        for pp in range(intgauss):
            detJ = Quad4.detJacobi(pt[pp], elementcoord)               
            B = PlaneStressIsoQuad4.getB(Model, elementcoord, pt[pp], nodedof) #np.dot(H, np.dot(invJ, diffN))
            BT = np.array(B).transpose() #transpose(B)
            BTC = np.dot(BT, C)
            K_elem_mat += np.dot(BTC, B)*L*np.abs(detJ)*wt[pp]*wt[pp] #multi_dot([BT, C, B])*L*detJ*wt[pp]*wt[pp] #dot(BTC, B)*L*detJ*wt[pp]*wt[pp]     
            # K_elem_mat += (np.transpose(B) @ C @ B)*L*detJ*wt[pp]*wt[pp]
        return K_elem_mat
    
    @cython.exceptval(check=False)
    @cython.boundscheck(False) # turn off bounds-checking for entire function
    @cython.wraparound(False)  # turn off negative index wrapping for entire function     
    def getMassConsistentMat(Model, double [:,:] inci, double [:,:] coord, double [:,:] tabmat, double [:,:] tabgeo, int intgauss, int element_number):
        cdef dict elem_set = PlaneStressIsoQuad4.getElementSet()
        cdef unsigned int nodedof = len(elem_set["dofs"]['d'])
        cdef dict shape_set = Quad4.getShapeSet()
        cdef unsigned int nodecon = len(shape_set['nodes'])
        cdef str type_shape = shape_set["key"]       
        cdef unsigned int edof = nodecon * nodedof
        cdef list nodelist = Quad4.getNodeList(inci, element_number)    
        cdef double [:,:] elementcoord = Quad4.getNodeCoord(coord, nodelist)
        cdef float E = tabmat[int(inci[element_number, 2]) - 1, 0]  # material elasticity
        cdef float v = tabmat[int(inci[element_number, 2]) - 1, 1]  # material poisson ratio
        cdef double [:,:] C = PlaneStressIsotropic.getElasticTensor(E, v)
        cdef float L = tabgeo[int(inci[element_number, 3] - 1), 4]
        cdef float R = tabmat[int(inci[element_number, 2]) - 1, 6]  # material density
        cdef list pt
        cdef list wt
        pt, wt = gauss_points(type_shape, intgauss) 
        cdef double [:,:] M_elem_mat = np.zeros((edof, edof), dtype=FLT64)
        cdef float detJ
        cdef double [:,:] N
        cdef double [:,:] NT
        cdef double [:,:] NTR
        for pp in range(intgauss):
            detJ = Quad4.detJacobi(pt[pp], elementcoord)
            N = Quad4.getShapeFunctions(pt[pp], nodedof)
            NT = N.transpose() #transpose(B)
            NTR = np.dot(NT, R)
            M_elem_mat += np.dot(NTR, N)*L*np.abs(detJ)*wt[pp]*wt[pp]  
        return M_elem_mat
    
    @cython.exceptval(check=False)
    @cython.boundscheck(False) # turn off bounds-checking for entire function
    @cython.wraparound(False)  # turn off negative index wrapping for entire function     
    def getElementDeformation(double [:] U, dict modelinfo):
        cdef int nodetot = modelinfo['nnode']
        cdef int nodedof = modelinfo['nodedof']
        cdef double [:,:] Udef = np.zeros((nodetot, 3), dtype=FLT64)
        cdef double [:,:] Umag = np.zeros((nodetot, 1), dtype=FLT64)
        cdef double [:,:] U_def = np.zeros((nodetot, 4), dtype=FLT64)
        for nn in range(1, nodetot + 1):
            Udef[nn - 1, 0] = U[nodedof * nn - 2]
            Udef[nn - 1, 1] = U[nodedof * nn - 1]
            Umag[nn - 1, 0] = np.sqrt(U[nodedof * nn - 2] ** 2 + U[nodedof * nn - 1] ** 2)
        U_def = np.concatenate((Umag, Udef), axis=1)
        return U_def


    def setTitleDeformation():
        return ["DISPL_X", "DISPL_Y", "DISPL_Z"] 
    
    @cython.exceptval(check=False)
    @cython.boundscheck(False) # turn off bounds-checking for entire function
    @cython.wraparound(False)  # turn off negative index wrapping for entire function     
    def getElementVolume(Model, double [:,:] inci, double [:,:] coord, double [:,:] tabgeo, int intgauss, int element_number):
        cdef float L = tabgeo[int(inci[element_number, 3] - 1), 4]
        cdef dict shape_set = Quad4.getShapeSet()
        cdef str type_shape = shape_set["key"]    
        cdef list nodelist = Quad4.getNodeList(inci, element_number)    
        cdef double [:,:] elementcoord = Quad4.getNodeCoord(coord, nodelist)
        cdef list pt
        cdef list wt
        pt, wt = gauss_points(type_shape, intgauss)
        cdef float detJ = 0.0
        cdef float Vol = 0.0
        for pp in range(intgauss):
            N = Quad4.N
            detJ += np.abs(Quad4.detJacobi(pt[pp], elementcoord))
        Vol = detJ*L
        return Vol