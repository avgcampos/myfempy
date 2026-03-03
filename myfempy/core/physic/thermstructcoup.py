from __future__ import annotations

import numpy as np
from scipy.special import roots_legendre

from myfempy.core.physic.loadstruct import LoadStructural
from myfempy.core.physic.structural import Structural
from myfempy.core.utilities import (gauss_points, get_elemen_from_nodelist,
                                    get_nodes_from_list, poly_area)


__docformat__ = "google"

__doc__ = """

==========================================================================
                            __                                
         _ __ ___   _   _  / _|  ___  _ __ ___   _ __   _   _ 
        | '_ ` _ \ | | | || |_  / _ \| '_ ` _ \ | '_ \ | | | |
        | | | | | || |_| ||  _||  __/| | | | | || |_) || |_| |
        |_| |_| |_| \__, ||_|   \___||_| |_| |_|| .__/  \__, |
                    |___/                       |_|     |___/ 
        myfempy -- MultiphYsics Finite Element Module to PYthon    
                    COMPUTATIONAL ANALYSIS PROGRAM                   
        Copyright (C) 2022-2026 Antonio Vinicius Garcia Campos        
==========================================================================
This Python file is part of myfempy project.

myfempy is a python package based on finite element method to multiphysics
analysis. The code is open source and *intended for educational and scientific
purposes only, not recommended to commercial use. The name myfempy is an acronym
for MultiphYsics Finite Elements Module to PYthon. You can help us by contributing
with the main project, send us a mensage on https://github.com/avgcampos/myfempy/discussions/10
If you use myfempy in your research, the  developers would be grateful if you 
could cite in your work.
																		
The code is written by Antonio Vinicius Garcia Campos.                                  
																		
A github repository, with the most up to date version of the code,      
can be found here: https://github.com/avgcampos/myfempy.                 
																		
The code is open source and intended for educational and scientific     
purposes only. If you use myfempy in your research, the developers      
would be grateful if you could cite this. The myfempy project is published
under the GPLv3, see the myfempy LICENSE on
https://github.com/avgcampos/myfempy/blob/main/LICENSE.
																		
Disclaimer:                                                             
The authors reserve all rights but do not guarantee that the code is    
free from errors. Furthermore, the authors shall not be liable in any   
event caused by the use of the program.

"""


class ThermalStructuralCoupling(Structural):
    """Thermal Structural Coupled field analysis <ConcreteClassService>"""

    def getLoadApply(Model, coupling):
        forcenodeaply = np.zeros((1, 4))
        if coupling["TYPE"] == "thermalstress":  # thermo stress mechanical
            fapp = ThermalStructuralCoupling.ForceThermalStress(Model, coupling)
            forcenodeaply = np.append(forcenodeaply, fapp, axis=0)
        else:
            pass
        forcenodeaply = forcenodeaply[1::][::]
        return forcenodeaply

    def ForceThermalStress(Model, coupling):
        forcenodedof = np.zeros((1, 4))
        inci = Model.inci
        coord = Model.coord
        tabmat = Model.tabmat
        tabgeo = Model.tabgeo
        intgauss = Model.intgauss
        
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
                fc_type_dof = np.tile([Model.modelinfo["dofs"]["f"]["fx"], Model.modelinfo["dofs"]["f"]["fy"], Model.modelinfo["dofs"]["f"]["fz"]], len(nodelist),)
                nodelist = np.repeat(nodelist, 3)
            except:
                fc_type_dof = np.tile([Model.modelinfo["dofs"]["f"]["fx"], Model.modelinfo["dofs"]["f"]["fy"]], len(nodelist),)
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
