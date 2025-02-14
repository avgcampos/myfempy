from __future__ import annotations

import numpy as np
from scipy.special import roots_legendre

from myfempy.core.physic.loadstruct import LoadStructural
from myfempy.core.physic.structural import Structural
from myfempy.core.utilities import (gauss_points, get_elemen_from_nodelist,
                                    get_nodes_from_list, poly_area)


class ThermalStructuralCoupling(Structural):
    """Thermal Structural Coupled field analysis <ConcreteClassService>"""

    def getLoadApply(Model, modelinfo, coupling):
        forcenodeaply = np.zeros((1, 4))
        if coupling["TYPE"] == "thermalstress":  # thermo stress mechanical
            fapp = ThermalStructuralCoupling.ForceThermalStress(
                Model, modelinfo, coupling
            )
            forcenodeaply = np.append(forcenodeaply, fapp, axis=0)
        else:
            pass
        forcenodeaply = forcenodeaply[1::][::]
        return forcenodeaply

    def ForceThermalStress(Model, modelinfo, coupling):
        forcenodedof = np.zeros((1, 4))
        inci = modelinfo["inci"]
        coord = modelinfo["coord"]
        tabmat = modelinfo["tabmat"]
        tabgeo = modelinfo["tabgeo"]
        intgauss = modelinfo["intgauss"]
        
        strain_thermal = Model.material.getStrainThermal(coupling["GRADTEMP"])

        for ee in range(inci.shape[0]):
            force_value_vector, nodelist = (
                ThermalStructuralCoupling.__body_thermal_stress(
                    Model,
                    inci,
                    coord,
                    tabmat,
                    tabgeo,
                    intgauss,
                    strain_thermal[:, ee],
                    ee,
                )
            )

            try:
                fc_type_dof = np.tile([modelinfo["dofs"]["f"]["fx"], modelinfo["dofs"]["f"]["fy"], modelinfo["dofs"]["f"]["fz"]], len(nodelist),)
                nodelist = np.repeat(nodelist, 3)
            except:
                fc_type_dof = np.tile([modelinfo["dofs"]["f"]["fx"], modelinfo["dofs"]["f"]["fy"]], len(nodelist),)
                nodelist = np.repeat(nodelist, 2)
                            
            for j in range(len(nodelist)):
                fcdof = np.array(
                    [
                        [
                            int(nodelist[j]),
                            fc_type_dof[j],
                            force_value_vector[j],
                            int(coupling["STEP"]),
                        ]
                    ]
                )
                forcenodedof = np.append(forcenodedof, fcdof, axis=0)
        # forcenodedof[np.nonzero(forcenodedof)]
        forcenodedof = forcenodedof[1::][::]
        return forcenodedof

    def getUpdateMatrix(Model, matrix, loadaply):
        return LoadStructural.getUpdateMatrix(Model, matrix, loadaply)

    def getUpdateLoad(self):
        return LoadStructural.getUpdateLoad(self)

    def __body_thermal_stress(
        Model,
        inci,
        coord,
        tabmat,
        tabgeo,
        intgauss,
        strain_thermal,
        element_number,
    ):
        # body force
        elem_set = Model.element.getElementSet()
        nodedof = len(elem_set["dofs"]["d"])
        shape_set = Model.shape.getShapeSet()
        nodecon = len(shape_set["nodes"])
        type_shape = shape_set["key"]
        edof = nodecon * nodedof
        nodelist = Model.shape.getNodeList(inci, element_number)
        elementcoord = Model.shape.getNodeCoord(coord, nodelist)
        a = tabmat[int(inci[element_number, 2]) - 1]["CTE"]
        C = Model.material.getElasticTensor(tabmat, inci, element_number)
        pt, wt = gauss_points(type_shape, intgauss)
        W = np.zeros((nodedof, 1))
        force_value_vector = np.zeros((edof, 1))
        for ip in range(intgauss):
            for jp in range(intgauss):
                for kp in range(intgauss):
                    detJ = Model.shape.getdetJacobi(np.array([pt[ip], pt[jp], pt[kp]]), elementcoord)
                    diffN = Model.shape.getDiffShapeFuntion(np.array([pt[ip], pt[jp], pt[kp]]), nodedof)
                    invJ = Model.shape.getinvJacobi(np.array([pt[ip], pt[jp], pt[kp]]), elementcoord, nodedof)
                    B = Model.element.getB(diffN, invJ)
                    force_value_vector += np.dot(np.dot(B.transpose(), C), strain_thermal.reshape((-1,1))) * (a * abs(detJ) * wt[ip] * wt[jp] * wt[kp])
        force_value_vector = np.reshape(force_value_vector, (edof))
        return force_value_vector, nodelist
