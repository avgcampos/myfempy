from __future__ import annotations

import numpy as np
from scipy.special import roots_legendre

from myfempy.core.physic.structural import Structural
from myfempy.core.utilities import (gauss_points, get_elemen_from_nodelist,
                                    get_nodes_from_list, poly_area)


class LoadStructural(Structural):
    """Structural Load Class <ConcreteClassService>"""

    def getLoadApply(Model, modelinfo, forcelist):
        forcenodeaply = np.zeros((1, 4))

        for fc_index in range(len(forcelist)):
            if forcelist[fc_index][0] == "forcenode":
                flist = forcelist[fc_index]
                fapp = LoadStructural.__ForceNodeLoadApply(modelinfo, flist)
                forcenodeaply = np.append(forcenodeaply, fapp, axis=0)

            elif forcelist[fc_index][0] == "forceedge":
                flist = forcelist[fc_index]
                fapp = LoadStructural.__ForceEdgeLoadApply(Model, modelinfo, flist)
                forcenodeaply = np.append(forcenodeaply, fapp, axis=0)

            elif forcelist[fc_index][0] == "forcesurf":
                flist = forcelist[fc_index]
                fapp = LoadStructural.__ForceSurfLoadApply(modelinfo, flist)
                forcenodeaply = np.append(forcenodeaply, fapp, axis=0)

            elif forcelist[fc_index][0] == "forcebody":
                flist = forcelist[fc_index]
                fapp = LoadStructural.__ForceBodyLoadApply(Model, modelinfo, flist)
                forcenodeaply = np.append(forcenodeaply, fapp, axis=0)
            else:
                pass

        forcenodeaply = forcenodeaply[1::][::]
        return forcenodeaply

    def setLoadDof(forcelist, node_list_fc):
        return Structural.setLoadDof(forcelist, node_list_fc)

    def __ForceNodeLoadApply(modelinfo, forcelist):
        forcenodedof = np.zeros((1, 4))

        nodelist = forcelist[3:]
        node_list_fc, dir_fc = get_nodes_from_list(
            nodelist, modelinfo["coord"], modelinfo["regions"]
        )

        force_value_vector = np.ones_like(node_list_fc) * float(forcelist[2])
        # fc_dof = forcelist[1]
        fc_type_dof = modelinfo["dofs"]["f"][forcelist[1]] * np.ones_like(node_list_fc)
        # fc_type_dof = LoadStructural.setLoadDof(fc_dof, node_list_fc)

        for j in range(len(node_list_fc)):
            fcdof = np.array(
                [
                    [
                        int(node_list_fc[j]),
                        fc_type_dof[j],
                        force_value_vector[j],
                        int(forcelist[8]),
                    ]
                ]
            )
            forcenodedof = np.append(forcenodedof, fcdof, axis=0)

        forcenodedof = forcenodedof[1::][::]
        return forcenodedof

    def __ForceBodyLoadApply(Model, modelinfo, forcelist):
        forcenodedof = np.zeros((1, 4))

        gravity_value = float(forcelist[2])

        inci = modelinfo["inci"]
        coord = modelinfo["coord"]
        tabmat = modelinfo["tabmat"]
        tabgeo = modelinfo["tabgeo"]
        intgauss = modelinfo["intgauss"]

        fc_type_dof = modelinfo["dofs"]["f"][forcelist[1]]

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
                            int(forcelist[8]),
                        ]
                    ]
                )
                forcenodedof = np.append(forcenodedof, fcdof, axis=0)

        forcenodedof = forcenodedof[1::][::]
        return forcenodedof

    def __ForceEdgeLoadApply(Model, modelinfo, forcelist):
        nodelist = forcelist[3:]
        node_list_fc, dir_fc = get_nodes_from_list(
            nodelist, modelinfo["coord"], modelinfo["regions"]
        )

        force_value = float(forcelist[2])
        force_dirc = forcelist[1]

        inci = modelinfo["inci"]
        coord = modelinfo["coord"]
        tabmat = modelinfo["tabmat"]
        tabgeo = modelinfo["tabgeo"]
        intgauss = modelinfo["intgauss"]

        fc_type = forcelist[1]

        elmlist = get_elemen_from_nodelist(inci, node_list_fc)

        forcenodedof = np.zeros((1, 4))
        for ee in range(len(elmlist)):
            (
                force_value_vector,
                nodelist,
                norm,
            ) = LoadStructural.__line_force_distribuition(
                Model,
                inci,
                coord,
                tabgeo,
                node_list_fc,
                elmlist[ee],
                force_value,
                fc_type,
            )

            if len(force_value_vector) > len(nodelist):
                nodelist = np.repeat(nodelist, len(modelinfo["dofs"]["f"]))
                fc_type_dof = np.tile(
                    [modelinfo["dofs"]["f"]["fx"], modelinfo["dofs"]["f"]["fy"]],
                    len(modelinfo["dofs"]["f"]),
                )

            elif forcelist[1] == "pressure":
                if int(norm[0]) == 1 and int(norm[1]) == 0:
                    fc_type_dof = modelinfo["dofs"]["f"]["fx"] * np.ones_like(nodelist)

                elif int(norm[0]) == 0 and int(norm[1]) == 1:
                    fc_type_dof = modelinfo["dofs"]["f"]["fy"] * np.ones_like(nodelist)

                else:
                    fc_type_dof = modelinfo["dofs"]["f"]["fy"] * np.ones_like(nodelist)

                # if get_side == '0' or get_side == '2':
                #     fc_type_dof = modelinfo['dofs']['f']['fy']*np.ones_like(nodelist)
                # elif get_side == '1' or get_side == '3':
                #     fc_type_dof = modelinfo['dofs']['f']['fx']*np.ones_like(nodelist)
                # else:
                #     fc_type_dof = modelinfo['dofs']['f']['fy']*np.ones_like(nodelist)

                # shape_set = Model.shape.getShapeSet()
                # fc_type_dof = np.tile(shape_set['sidenorm'][get_side], len(nodelist))

                # shape_set['sidenorm'][get_side]*np.ones_like(nodelist)

            else:
                fc_type_dof = modelinfo["dofs"]["f"][forcelist[1]] * np.ones_like(
                    nodelist
                )

            for j in range(len(nodelist)):
                fcdof = np.array(
                    [
                        [
                            int(nodelist[j]),
                            fc_type_dof[j],
                            force_value_vector[j],
                            int(forcelist[8]),
                        ]
                    ]
                )
                forcenodedof = np.append(forcenodedof, fcdof, axis=0)

        forcenodedof = forcenodedof[1::][::]
        return forcenodedof

    # def __ForceSurfLoadApply(modelinfo, forcelist):
    #     forcenodedof = np.zeros((1, 4))

    #     nodelist = forcelist[3:]
    #     node_list_fc, dir_fc = get_nodes_from_list(nodelist, modelinfo['coord'], modelinfo['regions'])

    #     force_value = float(forcelist[2])
    #     force_dirc = forcelist[1]

    #     force_value_vector = LoadStructural.__plane_force_distribuition(modelinfo['inci'], modelinfo['coord'], node_list_fc, force_value, force_dirc, modelinfo['elemid'], modelinfo['nodedof'])

    #     node_list_fc = np.repeat(node_list_fc, modelinfo["nodedof"], axis=0)

    #     fc_type_dof = modelinfo['dofs']['f'][forcelist[1]]*np.ones_like(node_list_fc)

    #     for j in range(len(node_list_fc)):

    #         fcdof = np.array(
    #             [[
    #                 int(node_list_fc[j]),
    #                 fc_type_dof[j],
    #                 force_value_vector[j, 0],
    #                 int(forcelist[8]),
    #                 ]])
    #         forcenodedof = np.append(forcenodedof, fcdof, axis=0)

    #     forcenodedof = forcenodedof[1::][::]
    #     return forcenodedof

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
        R = tabmat[int(inci[element_number, 2]) - 1, 6]  # material density
        t = tabgeo[int(inci[element_number, 3] - 1), 4]
        pt, wt = gauss_points(type_shape, intgauss)
        G = gravity_value
        W = np.zeros((nodedof, 1))
        W[fc_type_dof - 1, 0] = R * G
        force_value_vector = np.zeros((edof, 1))
        for pp in range(intgauss):
            detJ = Model.shape.getdetJacobi(pt[pp], elementcoord)
            N = Model.shape.getShapeFunctions(pt[pp], nodedof)
            force_value_vector += np.dot(N.transpose(), W) * t * abs(detJ) * wt[pp]
        force_value_vector = force_value_vector[np.nonzero(force_value_vector)]
        return force_value_vector, nodelist

    def __line_force_distribuition(
        Model, inci, coord, tabgeo, node_list_fc, elem, force_value, fc_type
    ):
        # edge load

        # elmlist = get_elemen_from_nodelist(inci, node_list_fc)

        elem_set = Model.element.getElementSet()
        nodedof = len(elem_set["dofs"]["d"])
        shape_set = Model.shape.getShapeSet()
        nodecon = len(shape_set["nodes"])
        type_shape = shape_set["key"]
        edof = nodecon * nodedof

        nodelist = Model.shape.getNodeList(inci, elem - 1)
        elementcoord = Model.shape.getNodeCoord(coord, nodelist)
        t = tabgeo[int(inci[elem - 1, 3] - 1), 4]

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

                dx = elementcoord[noj, 0] - elementcoord[noi, 0]
                dy = elementcoord[noj, 1] - elementcoord[noi, 1]

                L = np.sqrt(dx**2 + dy**2)

                tx = (-dy / L) * force_value
                ty = (dx / L) * force_value

                T = np.array([[tx], [ty]])  # force_value

                norm[0] = abs(dy / L)
                norm[1] = abs(dx / L)
            else:
                T = np.array([[0.0], [0.0]])

            idx_conec = np.array2string(idx_conec)
            get_side = LoadStructural.__get_side_fcapp(idx_conec[1:-1])

            infoside = Model.shape.getIsoParaSide(get_side)

            r_vl = infoside[0]
            r_ax = infoside[1]

            points, wt = gauss_points(type_shape, 2)

            points[:, r_ax] = r_vl

            force_value_vector = np.zeros((edof, 1))

            for pp in range(2):
                N = Model.shape.getShapeFunctions(points[pp], nodedof)

                diffN = Model.shape.getDiffShapeFuntion(points[pp], nodedof)
                J = Model.shape.getJacobian(points[pp], elementcoord)

                # detJ_e = np.sqrt(J[1,0]**2 + J[1,1]**2)

                detJ_e = Model.shape.getEdgeLength(J, get_side)

                force_value_vector += (
                    np.dot(N.transpose(), T) * t * abs(detJ_e) * wt[pp]
                )
            force_value_vector = force_value_vector[np.nonzero(force_value_vector)]
        return force_value_vector, nodes, norm

    def __get_side_fcapp(set_side):
        side = {
            "0 1": "0",  # quad4
            "1 0": "0",  # quad4
            "1 2": "1",  # quad4
            "2 1": "1",  # quad4
            "2 3": "2",  # quad4
            "3 2": "2",  # quad4
            "3 0": "3",  # quad4
            "0 3": "3",  # quad4
            "2 0": "2",  # tria3
            "0 2": "2",  # tria3
        }
        return side[set_side]

    # antigo

    # def __line_force_distribuition(inci, coord, tabgeo, node_list_fc, force_value, force_dirc, elemid, nodedof):

    #     elmlist = get_elemen_from_nodelist(inci, node_list_fc)

    #     if elemid == 2231:
    #         force_value_vector = np.zeros((nodedof * len(node_list_fc), 1))
    #         for i in range(len(elmlist)):
    #             noi = int(inci[int(elmlist[i] - 1), 4])
    #             noj = int(inci[int(elmlist[i] - 1), 5])
    #             nok = int(inci[int(elmlist[i] - 1), 6])
    #             if np.any([noi == node_list_fc[:]]):
    #                 no1 = noi
    #                 if np.any([noj == node_list_fc[:]]):
    #                     no2 = noj
    #                 elif np.any([nok == node_list_fc[:]]):
    #                     no2 = nok
    #                 else:
    #                     continue
    #             elif np.any([noj == node_list_fc[:]]):
    #                 no1 = noj
    #                 if np.any([nok == node_list_fc[:]]):
    #                     no2 = nok
    #                 else:
    #                     continue
    #             else:
    #                 continue

    #             L = np.sqrt(
    #                 (coord[no1 - 1, 1] - coord[no2 - 1, 1]) ** 2
    #                 + (coord[no1 - 1, 2] - coord[no2 - 1, 2]) ** 2
    #                 + (coord[no1 - 1, 3] - coord[no2 - 1, 3]) ** 2
    #             )
    #             no1dof = np.where(no1 == node_list_fc)[0][0]
    #             no2dof = np.where(no2 == node_list_fc)[0][0]
    #             loc = np.array(
    #                 [
    #                     nodedof * no1dof,
    #                     nodedof * no1dof + 1,
    #                     nodedof * no2dof,
    #                     nodedof * no2dof + 1,
    #                 ]
    #             )
    #             tck = tabgeo[
    #                 int(inci[int(elmlist[i] - 1), 3] - 1), 4
    #             ]

    #             if force_dirc == "fx":
    #                 force_value_vector[loc, 0] += np.array(
    #                     [force_value * tck * L / 2, 0, force_value * tck * L / 2, 0]
    #                 )

    #             elif force_dirc == "fy":
    #                 force_value_vector[loc, 0] += np.array(
    #                     [0, force_value * tck * L / 2, 0, force_value * tck * L / 2]
    #                 )

    #             # elif force_dirc == "pressure":

    #             #     # L = np.sqrt((nojx - noix) ** 2 + (nojy - noiy) ** 2)
    #             #     # s = (nojy - noiy) / L
    #             #     # c = (nojx - noix) / L

    #             #     S = (modelinfo["coord"][no2 - 1, 2] - modelinfo["coord"][no1 - 1, 2])/L
    #             #     C = (modelinfo["coord"][no2 - 1, 1] - modelinfo["coord"][no1 - 1, 1])/L

    #             #     force_value_vector[loc, 0] += np.array(
    #             #         [0, force_value * tck * L / 2, 0, force_value * tck * L / 2]
    #             #     )
    #             #     fc_type_dof = np.array([1, 2])
    #             else:
    #                 pass

    #     elif elemid == 2241:
    #         force_value_vector = np.zeros((nodedof * len(node_list_fc), 1))
    #         for i in range(len(elmlist)):
    #             noi = int(inci[int(elmlist[i] - 1), 4])
    #             noj = int(inci[int(elmlist[i] - 1), 5])
    #             nok = int(inci[int(elmlist[i] - 1), 6])
    #             nol = int(inci[int(elmlist[i] - 1), 7])
    #             if np.any([noi == node_list_fc[:]]):
    #                 no1 = noi
    #                 if np.any([noj == node_list_fc[:]]):
    #                     no2 = noj
    #                 elif np.any([nol == node_list_fc[:]]):
    #                     no2 = nol
    #                 else:
    #                     continue
    #             elif np.any([noj == node_list_fc[:]]):
    #                 no1 = noj
    #                 if np.any([nok == node_list_fc[:]]):
    #                     no2 = nok
    #                 else:
    #                     continue
    #             elif np.any([nok == node_list_fc[:]]):
    #                 no1 = nok
    #                 if np.any([nol == node_list_fc[:]]):
    #                     no2 = nol
    #                 else:
    #                     continue
    #             else:
    #                 continue

    #             L = np.sqrt(
    #                 (coord[no1 - 1, 1] - coord[no2 - 1, 1]) ** 2
    #                 + (coord[no1 - 1, 2] - coord[no2 - 1, 2]) ** 2
    #                 + (coord[no1 - 1, 3] - coord[no2 - 1, 3]) ** 2
    #             )
    #             no1dof = np.where(no1 == node_list_fc)[0][0]
    #             no2dof = np.where(no2 == node_list_fc)[0][0]

    #             loc = np.array(
    #                 [
    #                     nodedof * no1dof,
    #                     nodedof * no1dof + 1,
    #                     nodedof * no2dof,
    #                     nodedof * no2dof + 1,
    #                 ]
    #             )
    #             tck = tabgeo[
    #                 int(inci[int(elmlist[i] - 1), 3] - 1), 4
    #             ]

    #             if force_dirc == "fx":
    #                 force_value_vector[loc, 0] += np.array(
    #                     [force_value * tck * L / 2, 0, force_value * tck * L / 2, 0]
    #                 )

    #             elif force_dirc == "fy":
    #                 force_value_vector[loc, 0] += np.array(
    #                     [0, force_value * tck * L / 2, 0, force_value * tck * L / 2]
    #                 )

    #             # elif force_dirc == "pressure":

    #             #     # L = np.sqrt((nojx - noix) ** 2 + (nojy - noiy) ** 2)
    #             #     # s = (nojy - noiy) / L
    #             #     # c = (nojx - noix) / L

    #             #     print(no1, no2)
    #             #     # fvvecloc = np.array([0, -(force_value * tck * L / 2), 0, -(force_value * tck * L / 2)])

    #             #     vec1 = np.array([modelinfo["coord"][no1 - 1, 1], modelinfo["coord"][no1 - 1, 2]])
    #             #     vec2 = np.array([modelinfo["coord"][no2 - 1, 1], modelinfo["coord"][no2 - 1, 2]])

    #             #     angle  = angle_two_vector(vec1, vec2)

    #             #     alpha = angle_two_vector(np.array([1,0]), vec1)
    #             #     # phi = angle_two_vector(np.array([1,0]), vec1)

    #             #     beta = np.pi/2 - angle - alpha

    #             #     print('angle: ',np.degrees(angle),'alpha: ', np.degrees(alpha), 'beta: ', np.degrees(beta))

    #             #     S = (modelinfo["coord"][no2 - 1, 2] - modelinfo["coord"][no1 - 1, 2])/L
    #             #     C = (modelinfo["coord"][no2 - 1, 1] - modelinfo["coord"][no1 - 1, 1])/L

    #             #     force_value_vector[loc, 0] += np.array([-(force_value * tck * L / 2)*S*np.sign(S), -(force_value * tck * L / 2)*C*np.sign(C), -(force_value * tck * L / 2)*S*np.sign(S), -(force_value * tck * L / 2)*C*np.sign(C)])

    #             #     # Ro1X = (modelinfo["coord"][no1 - 1, 1])
    #             #     # Ro1Y = (modelinfo["coord"][no1 - 1, 2])

    #             #     # Ro1 = np.sqrt((Ro1X)**2 + (Ro1Y)**2)

    #             #     # C_the1 = Ro1X/Ro1
    #             #     # S_the1 = Ro1Y/Ro1

    #             #     # angle_rad1 = np.arctan(Ro1Y/Ro1X)
    #             #     # angle_deg1 = np.degrees(angle_rad1)

    #             #     # Ro2X = (modelinfo["coord"][no2 - 1, 1])
    #             #     # Ro2Y = (modelinfo["coord"][no2 - 1, 2])

    #             #     # angle_rad2 = np.arctan(Ro2Y/Ro2X)
    #             #     # angle_deg2 = np.degrees(angle_rad2)

    #             #     # Ro2 = np.sqrt((Ro2X)**2 + (Ro2Y)**2)

    #             #     # C_the2 = Ro2X/Ro2
    #             #     # S_the2 = Ro2Y/Ro2

    #             #     # T = np.zeros((4, 4))
    #             #     # T[0, 0] = c
    #             #     # T[0, 1] = -s
    #             #     # T[1, 0] = S_the1
    #             #     # T[1, 1] = C_the1
    #             #     # T[2, 2] = C_the2
    #             #     # T[2, 3] = -S_the2
    #             #     # T[3, 2] = S_the2
    #             #     # T[3, 3] = C_the2

    #             #     # force_value_vector[loc, 0] += np.dot(np.transpose(T), fvvecloc)
    #             #     # print(force_value_vector)
    #             #     print('======================')
    #             #     fc_type_dof = np.array([1, 2])

    #     elif elemid == 2341:
    #         force_value_vector = np.zeros((nodedof * len(node_list_fc), 1))
    #         for i in range(len(elmlist)):
    #             noi = int(inci[int(elmlist[i] - 1), 4])
    #             noj = int(inci[int(elmlist[i] - 1), 5])
    #             nok = int(inci[int(elmlist[i] - 1), 6])
    #             nol = int(inci[int(elmlist[i] - 1), 7])
    #             if np.any([noi == node_list_fc[:]]):
    #                 no1 = noi
    #                 if np.any([noj == node_list_fc[:]]):
    #                     no2 = noj
    #                 elif np.any([nol == node_list_fc[:]]):
    #                     no2 = nol
    #                 else:
    #                     continue
    #             elif np.any([noj == node_list_fc[:]]):
    #                 no1 = noj
    #                 if np.any([nok == node_list_fc[:]]):
    #                     no2 = nok
    #                 else:
    #                     continue
    #             elif np.any([nok == node_list_fc[:]]):
    #                 no1 = nok
    #                 if np.any([nol == node_list_fc[:]]):
    #                     no2 = nol
    #                 else:
    #                     continue
    #             else:
    #                 continue

    #             L = np.sqrt(
    #                 (coord[no1 - 1, 1] - coord[no2 - 1, 1]) ** 2
    #                 + (coord[no1 - 1, 2] - coord[no2 - 1, 2]) ** 2
    #                 + (coord[no1 - 1, 3] - coord[no2 - 1, 3]) ** 2
    #             )
    #             no1dof = np.where(no1 == node_list_fc)[0][0]
    #             no2dof = np.where(no2 == node_list_fc)[0][0]

    #             loc = np.array(
    #                 [
    #                     nodedof * no1dof,
    #                     nodedof * no1dof + 1,
    #                     nodedof * no1dof + 2,
    #                     nodedof * no2dof,
    #                     nodedof * no2dof + 1,
    #                     nodedof * no2dof + 2,
    #                 ]
    #             )
    #             tck = tabgeo[
    #                 int(inci[int(elmlist[i] - 1), 3] - 1), 4
    #             ]

    #             if force_dirc == "fz":
    #                 force_value_vector[loc, 0] += np.array(
    #                     [force_value * tck * L / 2, 0, 0, force_value * tck * L / 2, 0, 0]
    #                 )

    #             else:
    #                 pass

    #             # print(np.sum(force_value_vector))

    #     return force_value_vector

    # def __plane_force_distribuition(inci, coord, node_list_fc, force_value, force_dirc, elemid, nodedof):
    #     """force in surface appl.

    #     Arguments:
    #         modelinfo:dict        -- F.E. model dict with full information needed
    #         force_value:float     -- force value
    #         force_dirc:str        -- force direction
    #         node_list_fc:list     -- list of node with force applied
    #         fc_set:str            -- force set direction

    #     Returns:
    #         force_value_vector:np.array  -- force vecto
    #         fc_type_dof:list             -- force list dofs
    #     """
    #     # elmlist = np.array([0], dtype=int)
    #     # for ii in range(len(node_list_fc)):
    #     #     elm2list = inci[(np.asarray(np.where(inci[:, 4:] == node_list_fc[ii])))[0][:], 0,]
    #     #     elmlist = np.append(elmlist, elm2list)
    #     # elmlist = np.unique(elmlist)
    #     # elmlist = elmlist[1::][::]

    #     elmlist = get_elemen_from_nodelist(inci, node_list_fc)

    #     if elemid == 3381:
    #         force_value_vector = np.zeros((nodedof * len(node_list_fc), 1))
    #         for i in range(len(elmlist)):
    #             noi = int(inci[int(elmlist[i] - 1), 4])
    #             noj = int(inci[int(elmlist[i] - 1), 5])
    #             nok = int(inci[int(elmlist[i] - 1), 6])
    #             nol = int(inci[int(elmlist[i] - 1), 7])
    #             nom = int(inci[int(elmlist[i] - 1), 8])
    #             non = int(inci[int(elmlist[i] - 1), 9])
    #             noo = int(inci[int(elmlist[i] - 1), 10])
    #             nop = int(inci[int(elmlist[i] - 1), 11])
    #             if np.any([noi == node_list_fc[:]]):
    #                 no1 = noi
    #                 if np.any([noj == node_list_fc[:]]):
    #                     no2 = noj
    #                     if np.any([nok == node_list_fc[:]]):
    #                         no3 = nok
    #                         no4 = nol
    #                     else:
    #                         continue
    #                 elif np.any([nom == node_list_fc[:]]):
    #                     no2 = nom
    #                     if np.any([nop == node_list_fc[:]]):
    #                         no3 = nop
    #                         no4 = nol
    #                     elif np.any([non == node_list_fc[:]]):
    #                         no3 = non
    #                         no4 = noj
    #                     else:
    #                         continue
    #             elif np.any([nom == node_list_fc[:]]):
    #                 no1 = nom
    #                 if np.any([non == node_list_fc[:]]):
    #                     no2 = non
    #                     if np.any([noo == node_list_fc[:]]):
    #                         no3 = non
    #                         no4 = nop
    #                     else:
    #                         continue
    #                 else:
    #                     continue
    #             elif np.any([non == node_list_fc[:]]):
    #                 no1 = non
    #                 if np.any([noo == node_list_fc[:]]):
    #                     no2 = noo
    #                     no3 = nok
    #                     no4 = noj
    #                 else:
    #                     continue
    #             elif np.any([noo == node_list_fc[:]]):
    #                 no1 = noo
    #                 if np.any([nok == node_list_fc[:]]):
    #                     no2 = nok
    #                     no3 = nol
    #                     no4 = nop
    #                 else:
    #                     continue
    #             else:
    #                 continue
    #             coord_x_no1 = coord[no1 - 1, 1]
    #             coord_y_no1 = coord[no1 - 1, 2]
    #             coord_z_no1 = coord[no1 - 1, 3]
    #             coord_x_no2 = coord[no2 - 1, 1]
    #             coord_y_no2 = coord[no2 - 1, 2]
    #             coord_z_no2 = coord[no2 - 1, 3]
    #             coord_x_no3 = coord[no3 - 1, 1]
    #             coord_y_no3 = coord[no3 - 1, 2]
    #             coord_z_no3 = coord[no3 - 1, 3]
    #             coord_x_no4 = coord[no4 - 1, 1]
    #             coord_y_no4 = coord[no4 - 1, 2]
    #             coord_z_no4 = coord[no4 - 1, 3]
    #             poly = np.array(
    #                 [
    #                     [coord_x_no1, coord_y_no1, coord_z_no1],
    #                     [coord_x_no2, coord_y_no2, coord_z_no2],
    #                     [coord_x_no3, coord_y_no3, coord_z_no3],
    #                     [coord_x_no4, coord_y_no4, coord_z_no4],
    #                 ]
    #             )
    #             A = poly_area(poly)
    #             no1dof = np.where(no1 == node_list_fc)[0][0]
    #             no2dof = np.where(no2 == node_list_fc)[0][0]
    #             no3dof = np.where(no3 == node_list_fc)[0][0]
    #             no4dof = np.where(no4 == node_list_fc)[0][0]
    #             loc = np.array(
    #                 [
    #                     nodedof * no1dof,
    #                     nodedof * no1dof + 1,
    #                     nodedof * no1dof + 2,
    #                     nodedof * no2dof,
    #                     nodedof * no2dof + 1,
    #                     nodedof * no2dof + 2,
    #                     nodedof * no3dof,
    #                     nodedof * no3dof + 1,
    #                     nodedof * no3dof + 2,
    #                     nodedof * no4dof,
    #                     nodedof * no4dof + 1,
    #                     nodedof * no4dof + 2,
    #                 ]
    #             )
    #             if force_dirc == "fx":
    #                 force_value_vector[loc, 0] += np.array(
    #                     [
    #                         force_value * A / 4,
    #                         0.0,
    #                         0.0,
    #                         force_value * A / 4,
    #                         0.0,
    #                         0.0,
    #                         force_value * A / 4,
    #                         0.0,
    #                         0.0,
    #                         force_value * A / 4,
    #                         0.0,
    #                         0.0,
    #                     ]
    #                 )

    #             elif force_dirc == "fy":
    #                 force_value_vector[loc, 0] += np.array(
    #                     [
    #                         0.0,
    #                         force_value * A / 4,
    #                         0.0,
    #                         0.0,
    #                         force_value * A / 4,
    #                         0.0,
    #                         0.0,
    #                         force_value * A / 4,
    #                         0.0,
    #                         0.0,
    #                         force_value * A / 4,
    #                         0.0,
    #                     ]
    #                 )

    #             elif force_dirc == "fz":
    #                 force_value_vector[loc, 0] += np.array(
    #                     [
    #                         0.0,
    #                         0.0,
    #                         force_value * A / 4,
    #                         0.0,
    #                         0.0,
    #                         force_value * A / 4,
    #                         0.0,
    #                         0.0,
    #                         force_value * A / 4,
    #                         0.0,
    #                         0.0,
    #                         force_value * A / 4,
    #                     ]
    #                 )

    #     elif elemid == 3341:
    #         force_value_vector = np.zeros((nodedof * len(node_list_fc), 1))
    #         for i in range(len(elmlist)):
    #             noi = int(inci[int(elmlist[i] - 1), 4])
    #             noj = int(inci[int(elmlist[i] - 1), 5])
    #             nok = int(inci[int(elmlist[i] - 1), 6])
    #             nol = int(inci[int(elmlist[i] - 1), 7])
    #             if np.any([noi == node_list_fc[:]]):
    #                 no1 = noi
    #                 if np.any([noj == node_list_fc[:]]):
    #                     no2 = noj
    #                     if np.any([nok == node_list_fc[:]]):
    #                         no3 = nok
    #                     elif np.any([nol == node_list_fc[:]]):
    #                         no3 = nol
    #                     else:
    #                         continue
    #                 elif np.any([nok == node_list_fc[:]]):
    #                     no2 = nok
    #                     if np.any([nol == node_list_fc[:]]):
    #                         no3 = nol
    #                     else:
    #                         continue
    #                 else:
    #                     continue
    #             elif np.any([noj == node_list_fc[:]]):
    #                 no1 = noj
    #                 if np.any([nok == node_list_fc[:]]):
    #                     no2 = nok
    #                     if np.any([nol == node_list_fc[:]]):
    #                         no3 = nol
    #                     else:
    #                         continue
    #                 else:
    #                     continue
    #             else:
    #                 continue
    #             coord_x_no1 = coord[no1 - 1, 1]
    #             coord_y_no1 = coord[no1 - 1, 2]
    #             coord_z_no1 = coord[no1 - 1, 3]
    #             coord_x_no2 = coord[no2 - 1, 1]
    #             coord_y_no2 = coord[no2 - 1, 2]
    #             coord_z_no2 = coord[no2 - 1, 3]
    #             coord_x_no3 = coord[no3 - 1, 1]
    #             coord_y_no3 = coord[no3 - 1, 2]
    #             coord_z_no3 = coord[no3 - 1, 3]
    #             poly = np.array(
    #                 [
    #                     [coord_x_no1, coord_y_no1, coord_z_no1],
    #                     [coord_x_no2, coord_y_no2, coord_z_no2],
    #                     [coord_x_no3, coord_y_no3, coord_z_no3],
    #                 ]
    #             )
    #             A = poly_area(poly)
    #             no1dof = np.where(no1 == node_list_fc)[0][0]
    #             no2dof = np.where(no2 == node_list_fc)[0][0]
    #             no3dof = np.where(no3 == node_list_fc)[0][0]
    #             loc = np.array(
    #                 [
    #                     nodedof * no1dof,
    #                     nodedof * no1dof + 1,
    #                     nodedof * no1dof + 2,
    #                     nodedof * no2dof,
    #                     nodedof * no2dof + 1,
    #                     nodedof * no2dof + 2,
    #                     nodedof * no3dof,
    #                     nodedof * no3dof + 1,
    #                     nodedof * no3dof + 2,
    #                 ]
    #             )
    #             if force_dirc == "fx":
    #                 force_value_vector[loc, 0] += np.array(
    #                     [
    #                         force_value * A / 3,
    #                         0.0,
    #                         0.0,
    #                         force_value * A / 3,
    #                         0.0,
    #                         0.0,
    #                         force_value * A / 3,
    #                         0.0,
    #                         0.0,
    #                     ]
    #                 )
    #             elif force_dirc == "fy":
    #                 force_value_vector[loc, 0] += np.array(
    #                     [
    #                         0.0,
    #                         force_value * A / 3,
    #                         0.0,
    #                         0.0,
    #                         force_value * A / 3,
    #                         0.0,
    #                         0.0,
    #                         force_value * A / 3,
    #                         0.0,
    #                     ]
    #                 )
    #             elif force_dirc == "fz":
    #                 force_value_vector[loc, 0] += np.array(
    #                     [
    #                         0.0,
    #                         0.0,
    #                         force_value * A / 3,
    #                         0.0,
    #                         0.0,
    #                         force_value * A / 3,
    #                         0.0,
    #                         0.0,
    #                         force_value * A / 3,
    #                     ]
    #                 )
    #     return force_value_vector
