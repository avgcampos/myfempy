from __future__ import annotations

import numpy as np
from scipy.special import roots_legendre

from myfempy.core.physic.thermal import Thermal
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


class LoadThermal(Thermal):
    """Thermal Load Class <ConcreteClassService>"""

    def getLoadApply(Model, forcelist):
        forcenodeaply = np.zeros((1, 4))

        # if forcelist['TYPE'] == "forcenode":
        #     fapp = LoadThermal.__ForceNodeLoadApply(modelinfo, forcelist)
        #     forcenodeaply = np.append(forcenodeaply, fapp, axis=0)

        if forcelist["TYPE"] == "heatfluxedge":
            fapp = LoadThermal.ForceEdgeLoadApply(Model, forcelist)
            forcenodeaply = np.append(forcenodeaply, fapp, axis=0)

        elif forcelist['TYPE'] == "heatfluxsurf":
            fapp = LoadThermal.ForceSurfLoadApply(Model, forcelist)
            forcenodeaply = np.append(forcenodeaply, fapp, axis=0)

        elif forcelist["TYPE"] == "convectionedge":
            fapp = LoadThermal.ForceEdgeLoadApply(Model, forcelist)
            forcenodeaply = np.append(forcenodeaply, fapp, axis=0)

        elif forcelist['TYPE'] == "convectionsurf":
            fapp = LoadThermal.ForceSurfLoadApply(Model, forcelist)
            forcenodeaply = np.append(forcenodeaply, fapp, axis=0)

        elif forcelist["TYPE"] == "heatgeneration":
            fapp = LoadThermal.ForceBodyLoadApply(Model, forcelist)
            forcenodeaply = np.append(forcenodeaply, fapp, axis=0)

        else:
            pass

        forcenodeaply = forcenodeaply[1::][::]
        return forcenodeaply

    def ForceNodeLoadApply(Model, forcelist):
        forcenodedof = np.zeros((1, 4))
        nodelist = [
            forcelist["DIR"],
            forcelist["LOCX"],
            forcelist["LOCY"],
            forcelist["LOCZ"],
            forcelist["TAG"],
            forcelist["MESHNODE"],
        ]
        node_list_fc, dir_fc = get_nodes_from_list(
            nodelist, Model.coord,  Model.regions
        )
        force_value_vector = np.ones_like(node_list_fc) * float(forcelist["VAL"])
        fc_type_dof = Model.modelinfo["dofs"]["f"][forcelist["DOF"]] * np.ones_like(
            node_list_fc
        )
        for j in range(len(node_list_fc)):
            fcdof = np.array(
                [
                    [
                        int(node_list_fc[j]),
                        fc_type_dof[j],
                        force_value_vector[j],
                        int(forcelist["STEP"]),
                    ]
                ]
            )
            forcenodedof = np.append(forcenodedof, fcdof, axis=0)
        forcenodedof = forcenodedof[1::][::]
        return forcenodedof

    def ForceBodyLoadApply(Model, forcelist):
        forcenodedof = np.zeros((1, 4))
        heatgen = float(forcelist["VAL"])
        inci = Model.inci
        coord = Model.coord
        tabmat = Model.tabmat
        tabgeo = Model.tabgeo
        intgauss = Model.intgauss
        fc_type_dof = Model.modelinfo["dofs"]["f"][forcelist["DOF"]]
        for ee in range(inci.shape[0]):
            force_value_vector, nodelist = LoadThermal.__body_force_volumetric(
                Model,
                inci,
                coord,
                tabmat,
                tabgeo,
                intgauss,
                ee,
                heatgen,
                fc_type_dof,
            )
            for j in range(len(nodelist)):
                fcdof = np.array(
                    [
                        [
                            int(nodelist[j]),
                            fc_type_dof,
                            force_value_vector[j],
                            int(forcelist["STEP"]),
                        ]
                    ]
                )
                forcenodedof = np.append(forcenodedof, fcdof, axis=0)
        forcenodedof = forcenodedof[1::][::]
        return forcenodedof

    def ForceEdgeLoadApply(Model, forcelist):
        nodelist = [
            forcelist["DIR"],
            forcelist["LOCX"],
            forcelist["LOCY"],
            forcelist["LOCZ"],
            forcelist["TAG"],
            forcelist["MESHNODE"],
        ]  # forcelist[3:]
        node_list_fc, dir_fc = get_nodes_from_list(
            nodelist, Model.coord, Model.regions
        )
        force_value = float(forcelist["VAL"])
        force_dirc = forcelist["DOF"]
        inci = Model.inci
        coord = Model.coord
        tabmat = Model.tabmat
        tabgeo = Model.tabgeo
        intgauss = Model.intgauss
        fc_type = forcelist["DOF"]
        elmlist = get_elemen_from_nodelist(inci, node_list_fc)
        forcenodedof = np.zeros((1, 4))
        for ee in range(len(elmlist)):

            (
                force_value_vector,
                nodelist,
                norm,
            ) = LoadThermal.__line_force_distribuition(
                Model,
                inci,
                coord,
                tabgeo,
                intgauss,
                node_list_fc,
                elmlist[ee],
                force_value,
                fc_type,
            )
            
            fc_type_dof = Model.modelinfo["dofs"]["f"][forcelist["DOF"]] * np.ones_like(
                nodelist
            )

            for j in range(len(nodelist)):
                fcdof = np.array(
                    [
                        [
                            int(nodelist[j]),
                            fc_type_dof[j],
                            force_value_vector[j],
                            int(forcelist["STEP"]),
                        ]
                    ]
                )
                forcenodedof = np.append(forcenodedof, fcdof, axis=0)
        forcenodedof = forcenodedof[1::][::]
        return forcenodedof

    def ForceSurfLoadApply(Model, forcelist):
        nodelist = [
            forcelist["DIR"],
            forcelist["LOCX"],
            forcelist["LOCY"],
            forcelist["LOCZ"],
            forcelist["TAG"],
            forcelist["MESHNODE"],
        ]
        node_list_fc, dir_fc = get_nodes_from_list(
            nodelist, Model.coord, Model.regions
        )
        force_value = float(forcelist["VAL"])
        force_dirc = forcelist["DOF"]
        inci = Model.inci
        coord = Model.coord
        tabmat = Model.tabmat
        tabgeo = Model.tabgeo
        intgauss = Model.intgauss
        fc_type = forcelist["DOF"]
        elmlist = get_elemen_from_nodelist(inci, node_list_fc)
        forcenodedof = np.zeros((1, 4))
        for ee in range(len(elmlist)):
            (
                force_value_vector,
                nodelist,
                norm,
            ) = LoadThermal.__surf_force_distribuition(
                Model,
                inci,
                coord,
                tabgeo,
                intgauss,
                node_list_fc,
                elmlist[ee],
                force_value,
                fc_type,
            )
            
            fc_type_dof = Model.modelinfo["dofs"]["f"][forcelist["DOF"]] * np.ones_like(
                nodelist
            )

            for j in range(len(nodelist)):
                fcdof = np.array(
                    [
                        [
                            int(nodelist[j]),
                            fc_type_dof[j],
                            force_value_vector[j],
                            int(forcelist["STEP"]),
                        ]
                    ]
                )
                forcenodedof = np.append(forcenodedof, fcdof, axis=0)
        forcenodedof = forcenodedof[1::][::]
        return forcenodedof

    def getUpdateMatrix(Model, matrix, loadaply):
        addConv = np.where(loadaply[:, 1] == 15)
        if addConv[0].size:
            addLoad = loadaply[addConv, :][0]
            matrix["stiffness"] = Model.element.getUpdateMatrix(
                Model, matrix["stiffness"], addLoad
            )
        return matrix

    def getUpdateLoad(self):
        return None


    def __body_force_volumetric(
        Model,
        inci,
        coord,
        tabmat,
        tabgeo,
        intgauss,
        element_number,
        heat_gen,
        fc_type_dof,
    ):
        elem_set = Model.element.getElementSet()
        nodedof = len(elem_set["dofs"]["d"])
        shape_set = Model.shape.getShapeSet()
        nodecon = len(shape_set["nodes"])
        type_shape = shape_set["key"]
        edof = nodecon * nodedof
        nodelist = Model.shape.getNodeList(inci, element_number)
        elementcoord = Model.shape.getNodeCoord(coord, nodelist)
        pt, wt = gauss_points(type_shape, intgauss)
        Q = heat_gen
        force_value_vector = np.zeros((edof, 1))
        for ip in range(intgauss):
            for jp in range(intgauss):
                for kp in range(intgauss):
                    detJ = Model.shape.getdetJacobi(np.array([pt[ip], pt[jp], pt[kp]]), elementcoord)
                    N = Model.shape.getShapeFunctions(np.array([pt[ip], pt[jp], pt[kp]]), nodedof)
                    force_value_vector += np.dot(N.transpose(), Q) * abs(detJ) * wt[ip] * wt[jp] * wt[kp]
        force_value_vector = force_value_vector[np.nonzero(force_value_vector)]
        return force_value_vector, nodelist

    def __line_force_distribuition(
        Model,
        inci,
        coord,
        tabgeo,
        intgauss,
        node_list_fc,
        element_number,
        force_value,
        fc_type,
    ):
        elem_set = Model.element.getElementSet()
        nodedof = len(elem_set["dofs"]["d"])
        shape_set = Model.shape.getShapeSet()
        nodecon = len(shape_set["nodes"])
        type_shape = shape_set["key"]
        edof = nodecon * nodedof
        nodelist = Model.shape.getNodeList(inci, element_number - 1)
        elementcoord = Model.shape.getNodeCoord(coord, nodelist)
        t = tabgeo[int(inci[element_number - 1, 3] - 1)][
            "THICKN"
        ]  # tabgeo[int(inci[elem - 1, 3] - 1), 4]
        # nodes, idx_conec, __ = np.intersect1d(nodelist, node_list_fc, assume_unique=True, return_indices=True)
        test = np.in1d(nodelist, node_list_fc, assume_unique=True)
        nodes = np.array(nodelist)[test]
        idx_conec = np.where(test == True)[0]
        norm = np.zeros((2))
        q = force_value
        if len(idx_conec) < 2:
            nodes = np.repeat(nodes, 2)
            force_value_vector = np.zeros((len(nodes)))
            pass
        else:
            force_value_vector = np.zeros((len(nodes)))
            if force_value == 0:
                pass
            else:
                idx_conec = np.array2string(idx_conec)
                get_side = Model.shape.getSideAxis(idx_conec[1:-1])
                pt, wt = gauss_points(type_shape, intgauss)
                force_value_vector = np.zeros((edof, 1))
                for ip in range(intgauss):
                    points = Model.shape.getIsoParaSide(get_side, pt[ip])
                    N = Model.shape.getShapeFunctions(np.array(points), nodedof)
                    diffN = Model.shape.getDiffShapeFuntion(np.array(points), nodedof)
                    J = Model.shape.getJacobian(np.array(points), elementcoord)
                    detJ_e = Model.shape.getEdgeLength(J, get_side)
                    force_value_vector += (np.dot(N.transpose(), q) * t * abs(detJ_e) * wt[ip])
                force_value_vector = force_value_vector[np.nonzero(force_value_vector)]
        return force_value_vector, nodes, norm

    def __surf_force_distribuition(
        Model,
        inci,
        coord,
        tabgeo,
        intgauss,
        node_list_fc,
        element_number,
        force_value,
        fc_type,
    ):
        elem_set = Model.element.getElementSet()
        nodedof = len(elem_set["dofs"]["d"])
        shape_set = Model.shape.getShapeSet()
        nodecon = len(shape_set["nodes"])
        type_shape = shape_set["key"]
        edof = nodecon * nodedof
        nodelist = Model.shape.getNodeList(inci, element_number - 1)
        elementcoord = Model.shape.getNodeCoord(coord, nodelist)
        test = np.in1d(nodelist, node_list_fc, assume_unique=True)
        nodes = np.array(nodelist)[test]
        nodes_conec = np.where(test == True)[0]
        norm = np.zeros((3))
        q = force_value
        if len(nodes_conec) < 3:
            nodes = np.repeat(nodes, 2)
            force_value_vector = np.zeros((len(nodes)))
            pass
        else:
            force_value_vector = np.zeros((len(nodes)))
            if force_value == 0:
                pass
            else:
                idx_conec = np.array2string(nodes_conec)
                get_side = Model.shape.getSideAxis(idx_conec[1:-1])
                pt, wt = gauss_points(type_shape, intgauss)
                force_value_vector = np.zeros((edof, 1))
                for ip in range(intgauss):
                    for jp in range(intgauss):
                        points = Model.shape.getIsoParaSide(get_side, [pt[ip], pt[jp]]) # [pt[ip], pt[jp]]
                        N = Model.shape.getShapeFunctions(np.array(points), nodedof)
                        detJ_a = Model.shape.getAreaLength(get_side, elementcoord)
                        force_value_vector += (np.dot(N.transpose(), q) * abs(detJ_a) * wt[ip] * wt[jp])
                force_value_vector = force_value_vector[np.nonzero(force_value_vector)]
        return force_value_vector, nodes, norm