from __future__ import annotations

import numpy as np
from scipy.special import roots_legendre

from myfempy.core.physic.structural import Structural
from myfempy.core.utilities import (gauss_points, get_elemen_from_nodelist,
                                    get_nodes_from_list, poly_area, unit_normal)


class LoadStructural(Structural):
    """Structural Load Class <ConcreteClassService>"""

    def getLoadApply(Model, modelinfo, forcelist):
        forcenodeaply = np.zeros((1, 4))
        if forcelist['TYPE'] == "forcenode":
            fapp = LoadStructural.__ForceNodeLoadApply(modelinfo, forcelist)
            forcenodeaply = np.append(forcenodeaply, fapp, axis=0)
        elif forcelist['TYPE'] == "forceedge":
            fapp = LoadStructural.__ForceEdgeLoadApply(Model, modelinfo, forcelist)
            forcenodeaply = np.append(forcenodeaply, fapp, axis=0)
        elif forcelist['TYPE'] == "forcebeam":
            fapp = LoadStructural.__ForceBeamLoadApply(Model, modelinfo, forcelist)
            forcenodeaply = np.append(forcenodeaply, fapp, axis=0)
        elif forcelist['TYPE'] == "forcesurf":
            fapp = LoadStructural.__ForceSurfLoadApply(Model, modelinfo, forcelist)
            forcenodeaply = np.append(forcenodeaply, fapp, axis=0)
        elif forcelist['TYPE'] == "forcebody":
            fapp = LoadStructural.__ForceBodyLoadApply(Model, modelinfo, forcelist)
            forcenodeaply = np.append(forcenodeaply, fapp, axis=0)
        else:
            pass
        forcenodeaply = forcenodeaply[1::][::]
        return forcenodeaply

    def getUpdateMatrix(Model, matrix, loadaply):
        
        addSpring = np.where(loadaply[:, 1] == 16)
        addMass = np.where(loadaply[:, 1] == 15)
        
        if addSpring[0].size:
            addLoad = loadaply[addSpring, :][0]
            matrix["stiffness"] = Model.element.getUpdateMatrix(Model, matrix["stiffness"], addLoad)
       
        if addMass[0].size:
            addLoad = loadaply[addSpring, :][0]
            matrix["mass"] = Model.element.getUpdateMatrix(Model, matrix["mass"], addLoad)
        
        return matrix        

    def getUpdateLoad(self):
        return None

    def __ForceNodeLoadApply(modelinfo, forcelist):
        forcenodedof = np.zeros((1, 4))
        nodelist = [forcelist['DIR'], forcelist['LOCX'], forcelist['LOCY'], forcelist['LOCZ'], forcelist['TAG']]
        node_list_fc, dir_fc = get_nodes_from_list(
            nodelist, modelinfo["coord"], modelinfo["regions"]
        )
        force_value_vector = np.ones_like(node_list_fc) * float(forcelist['VAL'])
        fc_type_dof = modelinfo["dofs"]["f"][forcelist['DOF']] * np.ones_like(node_list_fc)
        for j in range(len(node_list_fc)):
            fcdof = np.array(
                [
                    [
                        int(node_list_fc[j]),
                        fc_type_dof[j],
                        force_value_vector[j],
                        int(forcelist['STEP']),
                    ]
                ]
            )
            forcenodedof = np.append(forcenodedof, fcdof, axis=0)
        forcenodedof = forcenodedof[1::][::]
        return forcenodedof

    def __ForceBodyLoadApply(Model, modelinfo, forcelist):
        forcenodedof = np.zeros((1, 4))
        gravity_value = float(forcelist['VAL'])
        inci = modelinfo["inci"]
        coord = modelinfo["coord"]
        tabmat = modelinfo["tabmat"]
        tabgeo = modelinfo["tabgeo"]
        intgauss = modelinfo["intgauss"]
        fc_type_dof = modelinfo["dofs"]["f"][forcelist['DOF']]
        for ee in range(inci.shape[0]):
            force_value_vector, nodelist = LoadStructural.__body_force_volumetric(
                Model,
                inci,
                coord,
                tabmat,
                tabgeo,
                intgauss,
                ee,
                gravity_value,
                fc_type_dof,
            )
            for j in range(len(nodelist)):
                fcdof = np.array(
                    [
                        [
                            int(nodelist[j]),
                            fc_type_dof,
                            force_value_vector[j],
                            int(forcelist['STEP']),
                        ]
                    ]
                )
                forcenodedof = np.append(forcenodedof, fcdof, axis=0)
        forcenodedof = forcenodedof[1::][::]
        return forcenodedof

    def __ForceEdgeLoadApply(Model, modelinfo, forcelist):
        nodelist = [forcelist['DIR'], forcelist['LOCX'], forcelist['LOCY'], forcelist['LOCZ'], forcelist['TAG']] #forcelist[3:]
        node_list_fc, dir_fc = get_nodes_from_list(
            nodelist, modelinfo["coord"], modelinfo["regions"]
        )
        force_value = float(forcelist['VAL'])
        force_dirc = forcelist['DOF']
        inci = modelinfo["inci"]
        coord = modelinfo["coord"]
        tabmat = modelinfo["tabmat"]
        tabgeo = modelinfo["tabgeo"]
        intgauss = modelinfo["intgauss"]
        fc_type = forcelist['DOF']
        elmlist = get_elemen_from_nodelist(inci, node_list_fc)
        forcenodedof = np.zeros((1, 4))
        for ee in range(len(elmlist)):
            force_value_vector, nodelist, norm = LoadStructural.__line_force_distribuition(
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
            
            elem_set = Model.element.getElementSet()
            nodedof = len(elem_set["dofs"]["d"])
                        
            if len(force_value_vector) > len(nodelist):
                nodelist = np.repeat(nodelist, nodedof)
                # fc_type_dof = np.tile([modelinfo["dofs"]["f"]["fx"], modelinfo["dofs"]["f"]["fy"]], int(len(nodelist)/nodedof))
                                
            if forcelist['DOF'] == "pressure":
                if int(norm[0]) == 1 and int(norm[1]) == 0:
                    fc_type_dof = modelinfo["dofs"]["f"]["fx"] * np.ones_like(nodelist)
                elif  int(norm[0]) == 0 and int(norm[1]) == 1:
                    fc_type_dof = modelinfo["dofs"]["f"]["fy"] * np.ones_like(nodelist)
                else:
                    fc_type_dof = np.tile([modelinfo["dofs"]["f"]["fx"], modelinfo["dofs"]["f"]["fy"]], int(len(nodelist)/nodedof)) #modelinfo["dofs"]["f"]["fy"] * np.ones_like(nodelist)
            else:
                fc_type_dof = modelinfo["dofs"]["f"][forcelist['DOF']] * np.ones_like(nodelist)
                
            for j in range(len(nodelist)):
                fcdof = np.array(
                    [
                        [
                            int(nodelist[j]),
                            fc_type_dof[j],
                            force_value_vector[j],
                            int(forcelist['STEP']),
                        ]
                    ]
                )
                forcenodedof = np.append(forcenodedof, fcdof, axis=0)
        forcenodedof = forcenodedof[1::][::]
        return forcenodedof

    def __ForceSurfLoadApply(Model, modelinfo, forcelist):
        nodelist = [forcelist['DIR'], forcelist['LOCX'], forcelist['LOCY'], forcelist['LOCZ'], forcelist['TAG']] #forcelist[3:]
        node_list_fc, dir_fc = get_nodes_from_list(
            nodelist, modelinfo["coord"], modelinfo["regions"]
        )
        force_value = float(forcelist['VAL'])
        force_dirc = forcelist['DOF']
        inci = modelinfo["inci"]
        coord = modelinfo["coord"]
        tabmat = modelinfo["tabmat"]
        tabgeo = modelinfo["tabgeo"]
        intgauss = modelinfo["intgauss"]
        fc_type = forcelist['DOF']
        elmlist = get_elemen_from_nodelist(inci, node_list_fc)
        forcenodedof = np.zeros((1, 4))
        for ee in range(len(elmlist)):
            force_value_vector, nodelist, norm = LoadStructural.__surf_force_distribuition(
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
            
            elem_set = Model.element.getElementSet()
            nodedof = len(elem_set["dofs"]["d"])
                        
            if len(force_value_vector) > len(nodelist):
                nodelist = np.repeat(nodelist, nodedof)
                # fc_type_dof = np.tile([modelinfo["dofs"]["f"]["fx"], modelinfo["dofs"]["f"]["fy"]], int(len(nodelist)/nodedof))
                                
            if forcelist['DOF'] == "pressure":
                if int(norm[0]) == 1 and int(norm[1]) == 0:
                    fc_type_dof = modelinfo["dofs"]["f"]["fx"] * np.ones_like(nodelist)
                elif  int(norm[0]) == 0 and int(norm[1]) == 1:
                    fc_type_dof = modelinfo["dofs"]["f"]["fy"] * np.ones_like(nodelist)
                else:
                    fc_type_dof = np.tile([modelinfo["dofs"]["f"]["fx"], modelinfo["dofs"]["f"]["fy"]], int(len(nodelist)/nodedof)) #modelinfo["dofs"]["f"]["fy"] * np.ones_like(nodelist)
            else:
                fc_type_dof = modelinfo["dofs"]["f"][forcelist['DOF']] * np.ones_like(nodelist)
                
            for j in range(len(nodelist)):
                fcdof = np.array(
                    [
                        [
                            int(nodelist[j]),
                            fc_type_dof[j],
                            force_value_vector[j],
                            int(forcelist['STEP']),
                        ]
                    ]
                )
                forcenodedof = np.append(forcenodedof, fcdof, axis=0)
        forcenodedof = forcenodedof[1::][::]
        return forcenodedof
    

    def __ForceBeamLoadApply(Model, modelinfo, forcelist):
        nodelist = [forcelist['DIR'], forcelist['LOCX'], forcelist['LOCY'], forcelist['LOCZ'], forcelist['TAG']] #forcelist[3:]
        node_list_fc, dir_fc = get_nodes_from_list(
            nodelist, modelinfo["coord"], modelinfo["regions"]
        )
        force_value = float(forcelist['VAL'])
        force_dirc = forcelist['DOF']
        inci = modelinfo["inci"]
        coord = modelinfo["coord"]
        tabmat = modelinfo["tabmat"]
        tabgeo = modelinfo["tabgeo"]
        intgauss = modelinfo["intgauss"]
        fc_type = forcelist['DOF']
        elmlist = get_elemen_from_nodelist(inci, node_list_fc)
        forcenodedof = np.zeros((1, 4))
        for ee in range(len(elmlist)):
            force_value_vector, nodeslist, norm = LoadStructural.__line_beam_distribuition(
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
            
            elem_set = Model.element.getElementSet()
            nodedof = len(elem_set["dofs"]["d"])
                        
            nodes = nodeslist
            if len(force_value_vector) > len(nodeslist):
                nodeslist = np.repeat(nodeslist, 2)
                # fc_type_dof = np.tile([modelinfo["dofs"]["f"]["fx"], modelinfo["dofs"]["f"]["fy"]], int(len(nodelist)/nodedof))
                                
            if forcelist['DOF'] == "pressure":
                if int(norm[0]) == 1 and int(norm[1]) == 0:
                    fc_type_dof = modelinfo["dofs"]["f"]["fx"] * np.ones_like(nodeslist)
                elif int(norm[0]) == 0 and int(norm[1]) == 1:
                    fc_type_dof = modelinfo["dofs"]["f"]["fy"] * np.ones_like(nodeslist)
                else:
                    fc_type_dof = np.tile([modelinfo["dofs"]["f"]["fx"], modelinfo["dofs"]["f"]["fy"]], int(len(nodes)/nodedof)) #modelinfo["dofs"]["f"]["fy"] * np.ones_like(nodelist)
            else:
                if force_dirc == "fy":
                    fc_type_dof = np.tile([2, 6], len(nodes))
                elif force_dirc == "fz":
                    fc_type_dof = np.tile([3, 5], len(nodes))
                elif force_dirc == "fx":
                    fc_type_dof = np.tile([1], len(nodes))
                elif force_dirc == "tx":
                    fc_type_dof = np.tile([4], len(nodes))
                else:
                    pass
                
            for j in range(len(nodeslist)):
                fcdof = np.array(
                    [
                        [
                            int(nodeslist[j]),
                            fc_type_dof[j],
                            force_value_vector[j],
                            int(forcelist['STEP']),
                        ]
                    ]
                )
                forcenodedof = np.append(forcenodedof, fcdof, axis=0)
        forcenodedof = forcenodedof[1::][::]
        return forcenodedof

    def __body_force_volumetric(
        Model,
        inci,
        coord,
        tabmat,
        tabgeo,
        intgauss,
        element_number,
        gravity_value,
        fc_type_dof,
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
        R = tabmat[int(inci[element_number, 2]) - 1]["RHO"] #tabmat[int(inci[element_number, 2]) - 1, 6]  # material density
        pt, wt = gauss_points(type_shape, intgauss)
        G = gravity_value
        W = np.zeros((nodedof, 1))
        W[fc_type_dof - 1, 0] = R * G
        force_value_vector = np.zeros((edof, 1))
        for ip in range(intgauss):
            for jp in range(intgauss):
                for kp in range(intgauss):
                    detJ = Model.shape.getdetJacobi(np.array([pt[ip], pt[jp], pt[kp]]), elementcoord)
                    N = Model.shape.getShapeFunctions(np.array([pt[ip], pt[jp], pt[kp]]), nodedof)
                    force_value_vector += np.dot(N.transpose(), W) * abs(detJ) * wt[ip] * wt[jp] * wt[kp]
        force_value_vector = force_value_vector[np.nonzero(force_value_vector)]
        # force_value_vector = np.reshape(force_value_vector, (edof))
        return force_value_vector, nodelist

    def __line_force_distribuition(
        Model, inci, coord, tabgeo, intgauss, node_list_fc, element_number, force_value, fc_type
    ):
        elem_set = Model.element.getElementSet()
        nodedof = len(elem_set["dofs"]["d"])
        shape_set = Model.shape.getShapeSet()
        nodecon = len(shape_set["nodes"])
        type_shape = shape_set["key"]
        edof = nodecon * nodedof
        nodelist = Model.shape.getNodeList(inci, element_number - 1)
        elementcoord = Model.shape.getNodeCoord(coord, nodelist)
        t = tabgeo[int(inci[element_number - 1, 3] - 1)]["THICKN"] #tabgeo[int(inci[elem - 1, 3] - 1), 4]
        # nodes, idx_conec, __ = np.intersect1d(nodelist, node_list_fc, assume_unique=True, return_indices=True)
        test = np.in1d(nodelist, node_list_fc, assume_unique=True)
        nodes = np.array(nodelist)[test]
        idx_conec = np.where(test == True)[0]
        norm = np.zeros((2))
        if len(idx_conec) < 2:
            nodes = np.repeat(nodes, 2)
            force_value_vector = np.zeros((len(nodes)))
            pass
        else:
            if fc_type == "fx":
                T = np.array([[force_value], [0.0]])  # force_value
                norm[0] = 1
            elif fc_type == "fy":
                T = np.array([[0.0], [force_value]])  # force_value
                norm[1] = 1
            elif fc_type == "pressure":  #  -->[+]<--
                noi = idx_conec[0]
                noj = idx_conec[1]
                dx = abs(elementcoord[noj, 0] - elementcoord[noi, 0])
                dy = abs(elementcoord[noj, 1] - elementcoord[noi, 1])
                L = np.sqrt(dx**2 + dy**2)
                # tx = (-dy / L) * force_value
                # ty = (dx / L) * force_value
                # T = np.array([[tx], [ty]])  # force_value
                norm[0] = dy / L
                norm[1] = dx / L
                T = -1*(np.array([norm]).T)*force_value #(np.array([norm]).T)*T
            else:
                T = np.array([[0.0], [0.0]])
            idx_conec = np.array2string(idx_conec)
            get_side = Model.shape.getSideAxis(idx_conec[1:-1]) #.__get_side_fcapp(idx_conec[1:-1])
            pt, wt = gauss_points(type_shape, intgauss)
            force_value_vector = np.zeros((edof, 1))
            for ip in range(intgauss):                
                points =  Model.shape.getIsoParaSide(get_side, pt[ip])
                N = Model.shape.getShapeFunctions(np.array(points), nodedof)
                diffN = Model.shape.getDiffShapeFuntion(np.array(points), nodedof)
                J = Model.shape.getJacobian(np.array(points), elementcoord)
                detJ_e = Model.shape.getEdgeLength(J, get_side)
                force_value_vector += np.dot(np.array(N).transpose(), T) * t * abs(detJ_e) * wt[ip]
            force_value_vector = force_value_vector[np.nonzero(force_value_vector)]
        return force_value_vector, nodes, norm
    
    def __surf_force_distribuition(
            Model, inci, coord, tabgeo, intgauss, node_list_fc, element_number, force_value, fc_type
        ):
            print(element_number)
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
            idx_conec = np.where(test == True)[0]
            norm = np.zeros((3))
            if len(idx_conec) < 3:
                nodes = np.repeat(nodes, 3)
                force_value_vector = np.zeros((len(nodes)))
                pass
            else:
                if fc_type == "fx":
                    T = np.array([[force_value], [0.0], [0.0]])  # force_value
                    norm[0] = 1
                elif fc_type == "fy":
                    T = np.array([[0.0], [force_value], [0.0]])  # force_value
                    norm[1] = 1
                elif fc_type == "fz":
                    T = np.array([[0.0], [0.0], [force_value]])  # force_value
                    norm[2] = 1
                # elif fc_type == "pressure":  #  -->[+]<--
                #     noi = idx_conec[0]
                #     noj = idx_conec[1]
                #     dx = abs(elementcoord[noj, 0] - elementcoord[noi, 0])
                #     dy = abs(elementcoord[noj, 1] - elementcoord[noi, 1])
                #     L = np.sqrt(dx**2 + dy**2)
                #     # tx = (-dy / L) * force_value
                #     # ty = (dx / L) * force_value
                #     # T = np.array([[tx], [ty]])  # force_value
                #     norm[0] = dy / L
                #     norm[1] = dx / L
                #     T = -1*(np.array([norm]).T)*force_value #(np.array([norm]).T)*T
                else:
                    T = np.array([[0.0], [0.0], [0.0]])
                idx_conec = np.array2string(idx_conec)
                get_side = Model.shape.getSideAxis(idx_conec[1:-1]) #.__get_side_fcapp(idx_conec[1:-1])
                pt, wt = gauss_points(type_shape, intgauss)
                force_value_vector = np.zeros((edof, 1))
                for ip in range(intgauss):
                    for jp in range(intgauss):                
                        points =  Model.shape.getIsoParaSide(get_side, [pt[ip], pt[jp]])
                        N = Model.shape.getShapeFunctions(np.array(points), nodedof)
                        diffN = Model.shape.getDiffShapeFuntion(np.array(points), nodedof)
                        J = Model.shape.getJacobian(np.array(points), elementcoord)
                        detJ_e = Model.shape.getAreaLength(J, get_side)
                        force_value_vector += np.dot(np.array(N).transpose(), T) * abs(detJ_e) * wt[ip] * wt[jp]
                force_value_vector = force_value_vector[np.nonzero(force_value_vector)]
            return force_value_vector, nodes, norm
    
    def __line_beam_distribuition(
        Model, inci, coord, tabgeo, intgauss, node_list_fc, element_number, force_value, fc_type
    ):
        elem_set = Model.element.getElementSet()
        nodedof = len(elem_set["dofs"]["d"])
        shape_set = Model.shape.getShapeSet()
        nodecon = len(shape_set["nodes"])
        type_shape = shape_set["key"]
        edof = nodecon * nodedof
        nodelist = Model.shape.getNodeList(inci, element_number - 1)
        elementcoord = Model.shape.getNodeCoord(coord, nodelist)
        t = tabgeo[int(inci[element_number - 1, 3] - 1)]["THICKN"] #tabgeo[int(inci[elem - 1, 3] - 1), 4]
        # nodes, idx_conec, __ = np.intersect1d(nodelist, node_list_fc, assume_unique=True, return_indices=True)
        test = np.in1d(nodelist, node_list_fc, assume_unique=True)
        nodes = np.array(nodelist)[test]
        idx_conec = np.where(test == True)[0]
        norm = np.zeros((3))
        if len(idx_conec) < 2:
            nodes = np.repeat(nodes, 2)
            force_value_vector = np.zeros((len(nodes)))
            pass
        else:
            if fc_type == "fx":
                T = np.array([[force_value], [0.0], [0.0], [0.0]])  # force_value
                norm[0] = 1
            elif fc_type == "fy":
                T = np.array([[0.0], [force_value], [0.0], [0.0]])  # force_value
                norm[1] = 1
            elif fc_type == "fz":
                T = np.array([[0.0], [0.0], [force_value], [0.0]])  # force_value
                norm[2] = 1
            elif fc_type == "tx":
                T = np.array([[0.0], [0.0], [0.0], [force_value]])  # force_value
                norm[0] = 1
            elif fc_type == "pressure":  #  -->[+]<--
                noi = idx_conec[0]
                noj = idx_conec[1]
                dx = abs(elementcoord[noj, 0] - elementcoord[noi, 0])
                dy = abs(elementcoord[noj, 1] - elementcoord[noi, 1])
                L = np.sqrt(dx**2 + dy**2)
                # tx = (-dy / L) * force_value
                # ty = (dx / L) * force_value
                # T = np.array([[tx], [ty]])  # force_value
                norm[0] = dy / L
                norm[1] = dx / L
                T = -1*(np.array([norm]).T)*force_value #(np.array([norm]).T)*T
            else:
                T = np.array([[0.0], [0.0]])
            idx_conec = np.array2string(idx_conec)
            get_side = Model.shape.getSideAxis(idx_conec[1:-1]) #.__get_side_fcapp(idx_conec[1:-1])
            pt, wt = gauss_points(type_shape, intgauss)
            force_value_vector = np.zeros((edof, 1))
            for ip in range(intgauss):                
                points =  Model.shape.getIsoParaSide(get_side, pt[ip])
                N = Model.shape.getShapeFunctions(np.array(points), nodedof)
                diffN = Model.shape.getDiffShapeFuntion(np.array(points), nodedof)
                J = Model.shape.getJacobian(np.array(points), elementcoord)
                detJ_e = Model.shape.getEdgeLength(J, get_side)
                force_value_vector += np.dot(np.array(N).transpose(), T) * abs(detJ_e) * wt[ip]
            force_value_vector = force_value_vector[np.nonzero(force_value_vector)]
        return force_value_vector, nodes, norm