import numpy as np

INT32 = np.uint32
FLT64 = np.float64

from myfempy.core.material.material import Material
from myfempy.core.utilities import get3D_LocalVector, getRotational_Matrix


class UniAxialStress(Material):
    """Uni-Axial Stress Isotropic Material Class <ConcreteClassService>"""

    def getMaterialSet():
        matset = {
            "mat": "uniaxialstress",
            "type": "isotropic",
        }
        return matset

    def getElasticTensor(Model=None, element_number=None):
        E = Model.tabmat[int(Model.inci[element_number, 2]) - 1]["EXX"]
        G = Model.tabmat[int(Model.inci[element_number, 2]) - 1]["GXY"]
        C = np.zeros((4, 4), dtype=FLT64)
        C[0, 0] = E
        C[1, 1] = E
        C[2, 2] = E
        C[3, 3] = G
        return C

    def getElementStrain(Model, U, ptg, element_number):
        elem_set = Model.element.getElementSet()
        nodedof = len(elem_set["dofs"]["d"])
        shape_set = Model.shape.getShapeSet()
        type_shape = shape_set["key"]

        nodelist = Model.shape.getNodeList(Model.inci, element_number)

        loc = Model.shape.getLocKey(nodelist, nodedof)

        elementcoord = Model.shape.getNodeCoord(Model.coord, nodelist)

        if type_shape == "line3":
            R = getRotational_Matrix(np.array(elementcoord[0:6]), 6)
        else:  # line2
            R = getRotational_Matrix(np.array(elementcoord[0:6]), 4)

        if type_shape == "line3":
            elementcoord_local = get3D_LocalVector(elementcoord, 3)
        else:
            elementcoord_local = get3D_LocalVector(elementcoord, 2)

        diffN = Model.shape.getDiffDiffShapeFuntion(np.array([ptg]), nodedof)
        
        invJ = Model.shape.getinvJacobi(np.array([ptg]), elementcoord, nodedof)
        
        B = Model.element.getB(diffN, invJ)

        cg = Model.geometry.getCGCoord(Model.tabgeo, Model.inci, element_number)

        epsilon = np.dot(B, np.dot(R, U[loc]))  # B @ (U[loc])

        strn_elm_normal_tension = epsilon[0]

        strn_elm_normal_bendingXY_max = -1.0 * cg["y_max"] * epsilon[1]
        
        strn_elm_normal_bendingXY_min = -1.0 * cg["y_min"] * epsilon[1]

        strn_elm_normal_bendingXZ_max = -1.0 * cg["z_max"] * epsilon[2]
        
        strn_elm_normal_bendingXZ_min = -1.0 * cg["z_min"] * epsilon[2]

        strn_elm_shear_torsion_max = cg["r_max"] * epsilon[3]

        strain = [
            strn_elm_normal_tension,
            strn_elm_normal_bendingXY_max,
            strn_elm_normal_bendingXY_min,
            strn_elm_normal_bendingXZ_max,
            strn_elm_normal_bendingXZ_min,
            strn_elm_shear_torsion_max,
            0.0,
        ]

        return epsilon, strain

    def getTitleStrain():
        title = [
            "STRAIN_NORMAL_TENSION",
            "STRAIN_NORMAL_BENDINGXY_MAX",
            "STRAIN_NORMAL_BENDINGXY_MIN",
            "STRAIN_NORMAL_BENDINGXZ_MAX",
            "STRAIN_NORMAL_BENDINGXZ_MIN",
            "STRAIN_SHEAR_TORSION",
            "null",
        ]
        return title

    def getElementStress(Model, epsilon, element_number):

        C = Model.material.getElasticTensor(Model, element_number)

        sigma = np.dot(C, epsilon)

        cg = Model.geometry.getCGCoord(Model.tabgeo, Model.inci, element_number)

        strs_elm_normal_tension = sigma[0]

        strs_elm_normal_bendingXY_max = -1.0 * cg["y_max"] * sigma[1]
        
        strs_elm_normal_bendingXY_min = -1.0 * cg["y_min"] * sigma[1]

        strs_elm_normal_bendingXZ_max = -1.0 * cg["z_max"] * sigma[2]
        
        strs_elm_normal_bendingXZ_min = -1.0 * cg["z_min"] * sigma[2]

        strs_elm_shear_torsion_max = cg["r_max"] * sigma[3]

        stress = [
            strs_elm_normal_tension,
            strs_elm_normal_bendingXY_max,
            strs_elm_normal_bendingXY_min,
            strs_elm_normal_bendingXZ_max,
            strs_elm_normal_bendingXZ_min,
            strs_elm_shear_torsion_max,
            0.0,
        ]

        return sigma, stress

    def getTitleStress():
        title = [
            "STRESS_NORMAL_TENSION",
            "STRESS_NORMAL_BENDINGXY_MAX",
            "STRESS_NORMAL_BENDINGXY_MIN",
            "STRESS_NORMAL_BENDINGXZ_MAX",
            "STRESS_NORMAL_BENDINGXZ_MIN",
            "STRESS_SHEAR_TORSION",
            "null",
        ]
        return title

    def getStrainEnergyDensity(sigma, epsilon, elemvol):
        strain_energy = 0.5 * np.dot(sigma.transpose(), epsilon)
        return strain_energy

    def getTitleCompliance():
        title = ["STRAIN_ENERGY_DENSITY"]
        return title

    def getFailureCriteria(sigma):
        return 0.0

    def getTitleFoS():
        title = ["FoS_YIELD_VON_MISES"]
        return title
