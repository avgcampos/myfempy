import numpy as np

INT32 = np.uint32
FLT64 = np.float64

from myfempy.core.material.material import Material


class HeatPlane(Material):
    """Heat Plane Isotropic Material Class <ConcreteClassService>"""

    def getMaterialSet():
        matset = {
            "mat": "heatplane",
            "type": "isotropic",
        }
        return matset

    def getElasticTensor(tabmat, inci, element_number, Model=None):
        Kxx = tabmat[int(inci[element_number, 2]) - 1]["KXX"]
        Kyy = tabmat[int(inci[element_number, 2]) - 1]["KYY"]
        D = np.zeros((2, 2), dtype=FLT64)
        D[0, 0] = Kxx
        D[1, 1] = Kyy
        return D

    def getElementGradTemp(Model, U, ptg, element_number):
        elem_set = Model.element.getElementSet()
        nodedof = len(elem_set["dofs"]["d"])

        nodelist = Model.shape.getNodeList(Model.inci, element_number)

        loc = Model.shape.getLocKey(nodelist, nodedof)

        elementcoord = Model.shape.getNodeCoord(Model.coord, nodelist)
        
        diffN = Model.shape.getDiffShapeFuntion(np.array([ptg, ptg]), nodedof)
        
        invJ = Model.shape.getinvJacobi(np.array([ptg, ptg]), elementcoord, nodedof)
        
        B = Model.element.getB(diffN, invJ)

        N = Model.shape.getShapeFunctions(np.array([ptg, ptg]), nodedof)

        epsilon = np.dot(B, U[loc])  # B @ (U[loc])

        epsilon_T = np.dot(N, U[loc])

        strn_elm_xx = epsilon[0]
        
        strn_elm_yy = epsilon[1]

        # strn_elm_vm = np.sqrt(epsilon[0]**2 + epsilon[1]**2)

        strain = [epsilon_T[0], strn_elm_xx, strn_elm_yy]

        return epsilon, strain

    def getTitleGradTemp():
        title = ["GRADTEMP", "GRADTEMP_XX", "GRADTEMP_YY"]
        return title

    def getElementHeatFlux(Model, epsilon, element_number):

        C = Model.material.getElasticTensor(Model.tabmat, Model.inci,  element_number)

        sigma = -1 * np.dot(C, epsilon)

        strs_elm_xx = sigma[0]
        
        strs_elm_yy = sigma[1]

        strs_elm_vm = np.sqrt(sigma[0] ** 2 + sigma[1] ** 2)

        stress = [strs_elm_vm, strs_elm_xx, strs_elm_yy]

        return sigma, stress

    def getTitleHeatFlux():
        title = ["HEATFLUX_MAG", "HEATFLUX_XX", "HEATFLUX_YY"]
        return title
