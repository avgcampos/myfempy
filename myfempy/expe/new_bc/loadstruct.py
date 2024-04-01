from __future__ import annotations

import numpy as np
# from scipy.special import roots_legendre

from myfempy.core.physic.structural import Structural
from myfempy.core.utilities import poly_area, get_nodes_from_list, get_elemen_from_nodelist
from myfempy.core.utilities import gauss_points

class LoadStructural(Structural):
    '''Structural Load Class <ConcreteClassService>'''

    def getLoadApply(Model, modelinfo, forcelist):
        forcenodeaply = np.zeros((1, 4))
                
        for fc_index in range(len(forcelist)):
            if forcelist[fc_index][0] == "forcenode":
                flist = forcelist[fc_index]
                fapp = LoadStructural.__ForceNodeLoadApply(modelinfo, flist)
                forcenodeaply = np.append(forcenodeaply, fapp, axis=0)
            
            elif forcelist[fc_index][0] == "forceedge":
                flist = forcelist[fc_index]
                fapp = LoadStructural.__ForceEdgeLoadApply(modelinfo, flist)
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
    
    def __ForceBodyLoadApply(Model, modelinfo, forcelist):
        forcenodedof = np.zeros((1, 4))      
       
        gravity_value = float(forcelist[2])

        inci = modelinfo['inci']
        coord = modelinfo['coord']
        tabmat = modelinfo['tabmat']
        tabgeo = modelinfo['tabgeo']
        intgauss = modelinfo['intgauss']

        fc_type_dof = modelinfo['dofs']['f'][forcelist[1]]

        for ee in range(inci.shape[0]):
            force_value_vector, nodelist = LoadStructural.__body_force_volumetric(Model, inci, coord, tabmat, tabgeo, intgauss, ee, gravity_value, fc_type_dof)
                        
            for j in range(len(nodelist)):
                fcdof = np.array(
                    [[
                        int(nodelist[j]),
                        fc_type_dof,
                        force_value_vector[j],
                        int(forcelist[7]),
                        ]])
                forcenodedof = np.append(forcenodedof, fcdof, axis=0)
        
        forcenodedof = forcenodedof[1::][::]
        return forcenodedof
    
    def __ForceNodeLoadApply(modelinfo, forcelist):

        forcenodedof = np.zeros((1, 4))      
       
        nodelist = forcelist[3:]
        node_list_fc, dir_fc = get_nodes_from_list(nodelist, modelinfo['coord'], modelinfo['regions'])

        force_value_vector = np.ones_like(node_list_fc) * float(forcelist[2])
        # fc_dof = forcelist[1]
        fc_type_dof = modelinfo['dofs']['f'][forcelist[1]]*np.ones_like(node_list_fc)
        # fc_type_dof = LoadStructural.setLoadDof(fc_dof, node_list_fc)

        for j in range(len(node_list_fc)):

            fcdof = np.array(
                [[
                    int(node_list_fc[j]),
                    fc_type_dof[j],
                    force_value_vector[j],
                    int(forcelist[7]),
                    ]])
            forcenodedof = np.append(forcenodedof, fcdof, axis=0)
        
        forcenodedof = forcenodedof[1::][::]
        return forcenodedof

    
    def __ForceEdgeLoadApply(modelinfo, forcelist):
        forcenodedof = np.zeros((1, 4))      
       
        nodelist = forcelist[3:]
        node_list_fc, dir_fc = get_nodes_from_list(nodelist, modelinfo['coord'], modelinfo['regions'])

        force_value = float(forcelist[2])
        force_dirc = forcelist[1]

        force_value_vector = LoadStructural.__line_force_distribuition(modelinfo['inci'], modelinfo['coord'], modelinfo['tabgeo'], node_list_fc, force_value, force_dirc, modelinfo['elemid'], modelinfo['nodedof'])
        

        node_list_fc = np.repeat(node_list_fc, modelinfo["nodedof"], axis=0)

        fc_type_dof = modelinfo['dofs']['f'][forcelist[1]]*np.ones_like(node_list_fc)       

        # fc_type_dof = LoadStructural.setLoadDof(forcelist[1], node_list_fc)
        for j in range(len(node_list_fc)):

            fcdof = np.array(
                [[
                    int(node_list_fc[j]),
                    fc_type_dof[j],
                    force_value_vector[j, 0],
                    int(forcelist[7]),
                    ]])
            forcenodedof = np.append(forcenodedof, fcdof, axis=0)
        
        forcenodedof = forcenodedof[1::][::]
        return forcenodedof
    

    def __ForceSurfLoadApply(modelinfo, forcelist):
        forcenodedof = np.zeros((1, 4))      
       
        nodelist = forcelist[3:]
        node_list_fc, dir_fc = get_nodes_from_list(nodelist, modelinfo['coord'], modelinfo['regions'])

        force_value = float(forcelist[2])
        force_dirc = forcelist[1]

        force_value_vector = LoadStructural.__plane_force_distribuition(modelinfo['inci'], modelinfo['coord'], node_list_fc, force_value, force_dirc, modelinfo['elemid'], modelinfo['nodedof'])
    
        node_list_fc = np.repeat(node_list_fc, modelinfo["nodedof"], axis=0)

        fc_type_dof = modelinfo['dofs']['f'][forcelist[1]]*np.ones_like(node_list_fc)       

        for j in range(len(node_list_fc)):

            fcdof = np.array(
                [[
                    int(node_list_fc[j]),
                    fc_type_dof[j],
                    force_value_vector[j, 0],
                    int(forcelist[7]),
                    ]])
            forcenodedof = np.append(forcenodedof, fcdof, axis=0)
        
        forcenodedof = forcenodedof[1::][::]
        return forcenodedof
      
      
    def __body_force_volumetric(Model, inci, coord, tabmat, tabgeo, intgauss, element_number, gravity_value, fc_type_dof):
        # body force
        elem_set = Model.element.getElementSet()
        nodedof = len(elem_set["dofs"]['d'])
        shape_set = Model.shape.getShapeSet()
        nodecon = len(shape_set['nodes'])
        type_shape = shape_set["key"]    
        edof = nodecon * nodedof
        nodelist = Model.shape.getNodeList(inci, element_number)    
        elementcoord = Model.shape.getNodeCoord(coord, nodelist)
        R = tabmat[int(inci[element_number, 2]) - 1, 6]  # material density
        t = tabgeo[int(inci[element_number, 3] - 1), 4]
        pt, wt = gauss_points(type_shape, intgauss)
        G = gravity_value
        W = np.zeros((nodedof,1))
        W[fc_type_dof-1, 0] = R*G
        force_value_vector = np.zeros((edof, 1))
        for pp in range(intgauss):
            detJ = Model.shape.getdetJacobi(pt[pp], elementcoord)
            N = Model.shape.getShapeFunctions(pt[pp], nodedof)            
            force_value_vector  += np.dot(N.transpose(), W)*t*abs(detJ)*wt[pp]
        force_value_vector = force_value_vector[np.nonzero(force_value_vector)]
        return force_value_vector, nodelist
      
    def __line_force_distribuition(inci, coord, tabgeo, node_list_fc, force_value, force_dirc, elemid, nodedof):
                
        elmlist = get_elemen_from_nodelist(inci, node_list_fc)

        if elemid == 2231:
            force_value_vector = np.zeros((nodedof * len(node_list_fc), 1))
            for i in range(len(elmlist)):
                noi = int(inci[int(elmlist[i] - 1), 4])
                noj = int(inci[int(elmlist[i] - 1), 5])
                nok = int(inci[int(elmlist[i] - 1), 6])
                if np.any([noi == node_list_fc[:]]):
                    no1 = noi
                    if np.any([noj == node_list_fc[:]]):
                        no2 = noj
                    elif np.any([nok == node_list_fc[:]]):
                        no2 = nok
                    else:
                        continue
                elif np.any([noj == node_list_fc[:]]):
                    no1 = noj
                    if np.any([nok == node_list_fc[:]]):
                        no2 = nok
                    else:
                        continue
                else:
                    continue
                               
                L = np.sqrt(
                    (coord[no1 - 1, 1] - coord[no2 - 1, 1]) ** 2
                    + (coord[no1 - 1, 2] - coord[no2 - 1, 2]) ** 2
                    + (coord[no1 - 1, 3] - coord[no2 - 1, 3]) ** 2
                )
                no1dof = np.where(no1 == node_list_fc)[0][0]
                no2dof = np.where(no2 == node_list_fc)[0][0]
                loc = np.array(
                    [
                        nodedof * no1dof,
                        nodedof * no1dof + 1,
                        nodedof * no2dof,
                        nodedof * no2dof + 1,
                    ]
                )
                tck = tabgeo[
                    int(inci[int(elmlist[i] - 1), 3] - 1), 4
                ]
                
                if force_dirc == "fx":
                    force_value_vector[loc, 0] += np.array(
                        [force_value * tck * L / 2, 0, force_value * tck * L / 2, 0]
                    )
                
                elif force_dirc == "fy":
                    force_value_vector[loc, 0] += np.array(
                        [0, force_value * tck * L / 2, 0, force_value * tck * L / 2]
                    )
                
                # elif force_dirc == "pressure":
                    
                #     # L = np.sqrt((nojx - noix) ** 2 + (nojy - noiy) ** 2)
                #     # s = (nojy - noiy) / L
                #     # c = (nojx - noix) / L
                    
                #     S = (modelinfo["coord"][no2 - 1, 2] - modelinfo["coord"][no1 - 1, 2])/L               
                #     C = (modelinfo["coord"][no2 - 1, 1] - modelinfo["coord"][no1 - 1, 1])/L    
                    
                #     force_value_vector[loc, 0] += np.array(
                #         [0, force_value * tck * L / 2, 0, force_value * tck * L / 2]
                #     )
                #     fc_type_dof = np.array([1, 2])
                else:
                    pass
                    
        elif elemid == 2241:
            force_value_vector = np.zeros((nodedof * len(node_list_fc), 1))
            for i in range(len(elmlist)):
                noi = int(inci[int(elmlist[i] - 1), 4])
                noj = int(inci[int(elmlist[i] - 1), 5])
                nok = int(inci[int(elmlist[i] - 1), 6])
                nol = int(inci[int(elmlist[i] - 1), 7])
                if np.any([noi == node_list_fc[:]]):
                    no1 = noi
                    if np.any([noj == node_list_fc[:]]):
                        no2 = noj
                    elif np.any([nol == node_list_fc[:]]):
                        no2 = nol
                    else:
                        continue
                elif np.any([noj == node_list_fc[:]]):
                    no1 = noj
                    if np.any([nok == node_list_fc[:]]):
                        no2 = nok
                    else:
                        continue
                elif np.any([nok == node_list_fc[:]]):
                    no1 = nok
                    if np.any([nol == node_list_fc[:]]):
                        no2 = nol
                    else:
                        continue
                else:
                    continue
                
                L = np.sqrt(
                    (coord[no1 - 1, 1] - coord[no2 - 1, 1]) ** 2
                    + (coord[no1 - 1, 2] - coord[no2 - 1, 2]) ** 2
                    + (coord[no1 - 1, 3] - coord[no2 - 1, 3]) ** 2
                )
                no1dof = np.where(no1 == node_list_fc)[0][0]
                no2dof = np.where(no2 == node_list_fc)[0][0]
                
                loc = np.array(
                    [
                        nodedof * no1dof,
                        nodedof * no1dof + 1,
                        nodedof * no2dof,
                        nodedof * no2dof + 1,
                    ]
                )
                tck = tabgeo[
                    int(inci[int(elmlist[i] - 1), 3] - 1), 4
                ]
                
                if force_dirc == "fx":
                    force_value_vector[loc, 0] += np.array(
                        [force_value * tck * L / 2, 0, force_value * tck * L / 2, 0]
                    )
                
                elif force_dirc == "fy":
                    force_value_vector[loc, 0] += np.array(
                        [0, force_value * tck * L / 2, 0, force_value * tck * L / 2]
                    )
                
                # elif force_dirc == "pressure":
                    
                #     # L = np.sqrt((nojx - noix) ** 2 + (nojy - noiy) ** 2)
                #     # s = (nojy - noiy) / L
                #     # c = (nojx - noix) / L
                    

                #     print(no1, no2)
                #     # fvvecloc = np.array([0, -(force_value * tck * L / 2), 0, -(force_value * tck * L / 2)])
                    
                #     vec1 = np.array([modelinfo["coord"][no1 - 1, 1], modelinfo["coord"][no1 - 1, 2]])
                #     vec2 = np.array([modelinfo["coord"][no2 - 1, 1], modelinfo["coord"][no2 - 1, 2]])

                #     angle  = angle_two_vector(vec1, vec2)

                #     alpha = angle_two_vector(np.array([1,0]), vec1)
                #     # phi = angle_two_vector(np.array([1,0]), vec1)

                #     beta = np.pi/2 - angle - alpha

                #     print('angle: ',np.degrees(angle),'alpha: ', np.degrees(alpha), 'beta: ', np.degrees(beta))

                #     S = (modelinfo["coord"][no2 - 1, 2] - modelinfo["coord"][no1 - 1, 2])/L              
                #     C = (modelinfo["coord"][no2 - 1, 1] - modelinfo["coord"][no1 - 1, 1])/L
                
                    
                #     force_value_vector[loc, 0] += np.array([-(force_value * tck * L / 2)*S*np.sign(S), -(force_value * tck * L / 2)*C*np.sign(C), -(force_value * tck * L / 2)*S*np.sign(S), -(force_value * tck * L / 2)*C*np.sign(C)])
                    
                #     # Ro1X = (modelinfo["coord"][no1 - 1, 1])
                #     # Ro1Y = (modelinfo["coord"][no1 - 1, 2])
                    
                #     # Ro1 = np.sqrt((Ro1X)**2 + (Ro1Y)**2)
                    
                #     # C_the1 = Ro1X/Ro1
                #     # S_the1 = Ro1Y/Ro1


                #     # angle_rad1 = np.arctan(Ro1Y/Ro1X)
                #     # angle_deg1 = np.degrees(angle_rad1)


                #     # Ro2X = (modelinfo["coord"][no2 - 1, 1])
                #     # Ro2Y = (modelinfo["coord"][no2 - 1, 2])

                #     # angle_rad2 = np.arctan(Ro2Y/Ro2X)
                #     # angle_deg2 = np.degrees(angle_rad2)
                    
                #     # Ro2 = np.sqrt((Ro2X)**2 + (Ro2Y)**2)
                    
                #     # C_the2 = Ro2X/Ro2
                #     # S_the2 = Ro2Y/Ro2


                #     # T = np.zeros((4, 4))
                #     # T[0, 0] = c
                #     # T[0, 1] = -s
                #     # T[1, 0] = S_the1
                #     # T[1, 1] = C_the1
                #     # T[2, 2] = C_the2
                #     # T[2, 3] = -S_the2
                #     # T[3, 2] = S_the2
                #     # T[3, 3] = C_the2
                                    
                #     # force_value_vector[loc, 0] += np.dot(np.transpose(T), fvvecloc)
                #     # print(force_value_vector)
                #     print('======================')
                #     fc_type_dof = np.array([1, 2])
        
        elif elemid == 2341:
            force_value_vector = np.zeros((nodedof * len(node_list_fc), 1))
            for i in range(len(elmlist)):
                noi = int(inci[int(elmlist[i] - 1), 4])
                noj = int(inci[int(elmlist[i] - 1), 5])
                nok = int(inci[int(elmlist[i] - 1), 6])
                nol = int(inci[int(elmlist[i] - 1), 7])
                if np.any([noi == node_list_fc[:]]):
                    no1 = noi
                    if np.any([noj == node_list_fc[:]]):
                        no2 = noj
                    elif np.any([nol == node_list_fc[:]]):
                        no2 = nol
                    else:
                        continue
                elif np.any([noj == node_list_fc[:]]):
                    no1 = noj
                    if np.any([nok == node_list_fc[:]]):
                        no2 = nok
                    else:
                        continue
                elif np.any([nok == node_list_fc[:]]):
                    no1 = nok
                    if np.any([nol == node_list_fc[:]]):
                        no2 = nol
                    else:
                        continue
                else:
                    continue
                
                L = np.sqrt(
                    (coord[no1 - 1, 1] - coord[no2 - 1, 1]) ** 2
                    + (coord[no1 - 1, 2] - coord[no2 - 1, 2]) ** 2
                    + (coord[no1 - 1, 3] - coord[no2 - 1, 3]) ** 2
                )
                no1dof = np.where(no1 == node_list_fc)[0][0]
                no2dof = np.where(no2 == node_list_fc)[0][0]
                
                loc = np.array(
                    [
                        nodedof * no1dof,
                        nodedof * no1dof + 1,
                        nodedof * no1dof + 2,
                        nodedof * no2dof,
                        nodedof * no2dof + 1,
                        nodedof * no2dof + 2,
                    ]
                )
                tck = tabgeo[
                    int(inci[int(elmlist[i] - 1), 3] - 1), 4
                ]
                
                if force_dirc == "fz":
                    force_value_vector[loc, 0] += np.array(
                        [force_value * tck * L / 2, 0, 0, force_value * tck * L / 2, 0, 0]
                    )
                
                else:
                    pass

                # print(np.sum(force_value_vector))

        return force_value_vector
    

    def __plane_force_distribuition(inci, coord, node_list_fc, force_value, force_dirc, elemid, nodedof):
        """force in surface appl.

        Arguments:
            modelinfo:dict        -- F.E. model dict with full information needed
            force_value:float     -- force value
            force_dirc:str        -- force direction
            node_list_fc:list     -- list of node with force applied
            fc_set:str            -- force set direction

        Returns:
            force_value_vector:np.array  -- force vecto
            fc_type_dof:list             -- force list dofs
        """
        # elmlist = np.array([0], dtype=int)
        # for ii in range(len(node_list_fc)):
        #     elm2list = inci[(np.asarray(np.where(inci[:, 4:] == node_list_fc[ii])))[0][:], 0,]
        #     elmlist = np.append(elmlist, elm2list)
        # elmlist = np.unique(elmlist)
        # elmlist = elmlist[1::][::]
        
        elmlist = get_elemen_from_nodelist(inci, node_list_fc)
       
        if elemid == 3381:
            force_value_vector = np.zeros((nodedof * len(node_list_fc), 1))
            for i in range(len(elmlist)):
                noi = int(inci[int(elmlist[i] - 1), 4])
                noj = int(inci[int(elmlist[i] - 1), 5])
                nok = int(inci[int(elmlist[i] - 1), 6])
                nol = int(inci[int(elmlist[i] - 1), 7])
                nom = int(inci[int(elmlist[i] - 1), 8])
                non = int(inci[int(elmlist[i] - 1), 9])
                noo = int(inci[int(elmlist[i] - 1), 10])
                nop = int(inci[int(elmlist[i] - 1), 11])
                if np.any([noi == node_list_fc[:]]):
                    no1 = noi
                    if np.any([noj == node_list_fc[:]]):
                        no2 = noj
                        if np.any([nok == node_list_fc[:]]):
                            no3 = nok
                            no4 = nol
                        else:
                            continue
                    elif np.any([nom == node_list_fc[:]]):
                        no2 = nom
                        if np.any([nop == node_list_fc[:]]):
                            no3 = nop
                            no4 = nol
                        elif np.any([non == node_list_fc[:]]):
                            no3 = non
                            no4 = noj
                        else:
                            continue
                elif np.any([nom == node_list_fc[:]]):
                    no1 = nom
                    if np.any([non == node_list_fc[:]]):
                        no2 = non
                        if np.any([noo == node_list_fc[:]]):
                            no3 = non
                            no4 = nop
                        else:
                            continue
                    else:
                        continue
                elif np.any([non == node_list_fc[:]]):
                    no1 = non
                    if np.any([noo == node_list_fc[:]]):
                        no2 = noo
                        no3 = nok
                        no4 = noj
                    else:
                        continue
                elif np.any([noo == node_list_fc[:]]):
                    no1 = noo
                    if np.any([nok == node_list_fc[:]]):
                        no2 = nok
                        no3 = nol
                        no4 = nop
                    else:
                        continue
                else:
                    continue
                coord_x_no1 = coord[no1 - 1, 1]
                coord_y_no1 = coord[no1 - 1, 2]
                coord_z_no1 = coord[no1 - 1, 3]
                coord_x_no2 = coord[no2 - 1, 1]
                coord_y_no2 = coord[no2 - 1, 2]
                coord_z_no2 = coord[no2 - 1, 3]
                coord_x_no3 = coord[no3 - 1, 1]
                coord_y_no3 = coord[no3 - 1, 2]
                coord_z_no3 = coord[no3 - 1, 3]
                coord_x_no4 = coord[no4 - 1, 1]
                coord_y_no4 = coord[no4 - 1, 2]
                coord_z_no4 = coord[no4 - 1, 3]
                poly = np.array(
                    [
                        [coord_x_no1, coord_y_no1, coord_z_no1],
                        [coord_x_no2, coord_y_no2, coord_z_no2],
                        [coord_x_no3, coord_y_no3, coord_z_no3],
                        [coord_x_no4, coord_y_no4, coord_z_no4],
                    ]
                )
                A = poly_area(poly)
                no1dof = np.where(no1 == node_list_fc)[0][0]
                no2dof = np.where(no2 == node_list_fc)[0][0]
                no3dof = np.where(no3 == node_list_fc)[0][0]
                no4dof = np.where(no4 == node_list_fc)[0][0]
                loc = np.array(
                    [
                        nodedof * no1dof,
                        nodedof * no1dof + 1,
                        nodedof * no1dof + 2,
                        nodedof * no2dof,
                        nodedof * no2dof + 1,
                        nodedof * no2dof + 2,
                        nodedof * no3dof,
                        nodedof * no3dof + 1,
                        nodedof * no3dof + 2,
                        nodedof * no4dof,
                        nodedof * no4dof + 1,
                        nodedof * no4dof + 2,
                    ]
                )
                if force_dirc == "fx":
                    force_value_vector[loc, 0] += np.array(
                        [
                            force_value * A / 4,
                            0.0,
                            0.0,
                            force_value * A / 4,
                            0.0,
                            0.0,
                            force_value * A / 4,
                            0.0,
                            0.0,
                            force_value * A / 4,
                            0.0,
                            0.0,
                        ]
                    )

                elif force_dirc == "fy":
                    force_value_vector[loc, 0] += np.array(
                        [
                            0.0,
                            force_value * A / 4,
                            0.0,
                            0.0,
                            force_value * A / 4,
                            0.0,
                            0.0,
                            force_value * A / 4,
                            0.0,
                            0.0,
                            force_value * A / 4,
                            0.0,
                        ]
                    )

                elif force_dirc == "fz":
                    force_value_vector[loc, 0] += np.array(
                        [
                            0.0,
                            0.0,
                            force_value * A / 4,
                            0.0,
                            0.0,
                            force_value * A / 4,
                            0.0,
                            0.0,
                            force_value * A / 4,
                            0.0,
                            0.0,
                            force_value * A / 4,
                        ]
                    )
                    
        elif elemid == 3341:
            force_value_vector = np.zeros((nodedof * len(node_list_fc), 1))
            for i in range(len(elmlist)):
                noi = int(inci[int(elmlist[i] - 1), 4])
                noj = int(inci[int(elmlist[i] - 1), 5])
                nok = int(inci[int(elmlist[i] - 1), 6])
                nol = int(inci[int(elmlist[i] - 1), 7])
                if np.any([noi == node_list_fc[:]]):
                    no1 = noi
                    if np.any([noj == node_list_fc[:]]):
                        no2 = noj
                        if np.any([nok == node_list_fc[:]]):
                            no3 = nok
                        elif np.any([nol == node_list_fc[:]]):
                            no3 = nol
                        else:
                            continue
                    elif np.any([nok == node_list_fc[:]]):
                        no2 = nok
                        if np.any([nol == node_list_fc[:]]):
                            no3 = nol
                        else:
                            continue
                    else:
                        continue
                elif np.any([noj == node_list_fc[:]]):
                    no1 = noj
                    if np.any([nok == node_list_fc[:]]):
                        no2 = nok
                        if np.any([nol == node_list_fc[:]]):
                            no3 = nol
                        else:
                            continue
                    else:
                        continue
                else:
                    continue
                coord_x_no1 = coord[no1 - 1, 1]
                coord_y_no1 = coord[no1 - 1, 2]
                coord_z_no1 = coord[no1 - 1, 3]
                coord_x_no2 = coord[no2 - 1, 1]
                coord_y_no2 = coord[no2 - 1, 2]
                coord_z_no2 = coord[no2 - 1, 3]
                coord_x_no3 = coord[no3 - 1, 1]
                coord_y_no3 = coord[no3 - 1, 2]
                coord_z_no3 = coord[no3 - 1, 3]
                poly = np.array(
                    [
                        [coord_x_no1, coord_y_no1, coord_z_no1],
                        [coord_x_no2, coord_y_no2, coord_z_no2],
                        [coord_x_no3, coord_y_no3, coord_z_no3],
                    ]
                )
                A = poly_area(poly)
                no1dof = np.where(no1 == node_list_fc)[0][0]
                no2dof = np.where(no2 == node_list_fc)[0][0]
                no3dof = np.where(no3 == node_list_fc)[0][0]
                loc = np.array(
                    [
                        nodedof * no1dof,
                        nodedof * no1dof + 1,
                        nodedof * no1dof + 2,
                        nodedof * no2dof,
                        nodedof * no2dof + 1,
                        nodedof * no2dof + 2,
                        nodedof * no3dof,
                        nodedof * no3dof + 1,
                        nodedof * no3dof + 2,
                    ]
                )
                if force_dirc == "fx":
                    force_value_vector[loc, 0] += np.array(
                        [
                            force_value * A / 3,
                            0.0,
                            0.0,
                            force_value * A / 3,
                            0.0,
                            0.0,
                            force_value * A / 3,
                            0.0,
                            0.0,
                        ]
                    )
                elif force_dirc == "fy":
                    force_value_vector[loc, 0] += np.array(
                        [
                            0.0,
                            force_value * A / 3,
                            0.0,
                            0.0,
                            force_value * A / 3,
                            0.0,
                            0.0,
                            force_value * A / 3,
                            0.0,
                        ]
                    )
                elif force_dirc == "fz":
                    force_value_vector[loc, 0] += np.array(
                        [
                            0.0,
                            0.0,
                            force_value * A / 3,
                            0.0,
                            0.0,
                            force_value * A / 3,
                            0.0,
                            0.0,
                            force_value * A / 3,
                        ]
                    )
        return force_value_vector