#!/usr/bin/env python
__doc__ ="""
forces to nodes
"""
import numpy as np


def force_edge(modelinfo, force_value, force_dirc, node_list_fc, dir_fc):
    elmlist = [None]
    for ii in range(len(node_list_fc)):
        elm2list = modelinfo['inci'][(np.asarray(
            np.where(modelinfo['inci'][:, 4:] == node_list_fc[ii])))[0][:], 0]
        elmlist.extend(elm2list)
    elmlist = elmlist[1::][::]
    elmlist = np.unique(elmlist)
    if (dir_fc == 'x') or (dir_fc == 'y_x') or (dir_fc == 'z_x'):
        coord_fc = 1
    elif (dir_fc == 'y') or (dir_fc == 'x_y') or (dir_fc == 'z_y'):
        coord_fc = 2
    elif (dir_fc == 'z') or (dir_fc == 'x_z') or (dir_fc == 'y_z'):
        coord_fc = 3
    if modelinfo['elemid'][0] == 210:
        force_value_vector = np.zeros(
            (modelinfo['nodedof'][0]*len(node_list_fc), 1))
        for i in range(len(elmlist)):
            noi = int(modelinfo['inci'][int(elmlist[i]-1), 4])
            noj = int(modelinfo['inci'][int(elmlist[i]-1), 5])
            nok = int(modelinfo['inci'][int(elmlist[i]-1), 6])
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
            L = np.sqrt((modelinfo['coord'][no1-1, 1] - modelinfo['coord'][no2-1, 1])**2 + (modelinfo['coord'][no1-1, 2] - modelinfo['coord'][no2-1, 2])**2 +
                        (modelinfo['coord'][no1-1, 3] - modelinfo['coord'][no2-1, 3])**2)
            no1dof = np.where(no1 == node_list_fc)[0][0]
            no2dof = np.where(no2 == node_list_fc)[0][0]
            loc = np.array([modelinfo['nodedof'][0]*no1dof, modelinfo['nodedof'][0] *
                           no1dof+1, modelinfo['nodedof'][0]*no2dof, modelinfo['nodedof'][0]*no2dof+1])
            tck = modelinfo['tabgeo'][int(
                modelinfo['inci'][int(elmlist[i]-1), 3]-1), 4]
            if force_dirc == 'fx':
                force_value_vector[loc, 0] += np.array(
                    [force_value*tck*L/2, 0, force_value*tck*L/2, 0])
                fc_type_dof = np.array([1, 0])
            elif force_dirc == 'fy':
                force_value_vector[loc, 0] += np.array(
                    [0, force_value*tck*L/2, 0, force_value*tck*L/2])
                fc_type_dof = np.array([0, 2])
    elif modelinfo['elemid'][0] == 220:
        force_value_vector = np.zeros(
            (modelinfo['nodedof'][0]*len(node_list_fc), 1))
        for i in range(len(elmlist)):
            noi = int(modelinfo['inci'][int(elmlist[i]-1), 4])
            noj = int(modelinfo['inci'][int(elmlist[i]-1), 5])
            nok = int(modelinfo['inci'][int(elmlist[i]-1), 6])
            nol = int(modelinfo['inci'][int(elmlist[i]-1), 7])
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
            L = np.sqrt((modelinfo['coord'][no1-1, 1] - modelinfo['coord'][no2-1, 1])**2 + (modelinfo['coord'][no1-1, 2] - modelinfo['coord'][no2-1, 2])**2 +
                        (modelinfo['coord'][no1-1, 3] - modelinfo['coord'][no2-1, 3])**2)
            no1dof = np.where(no1 == node_list_fc)[0][0]
            no2dof = np.where(no2 == node_list_fc)[0][0]
            loc = np.array([modelinfo['nodedof'][0]*no1dof, modelinfo['nodedof'][0] *
                           no1dof+1, modelinfo['nodedof'][0]*no2dof, modelinfo['nodedof'][0]*no2dof+1])
            tck = modelinfo['tabgeo'][int(
                modelinfo['inci'][int(elmlist[i]-1), 3]-1), 4]
            if force_dirc == 'fx':
                force_value_vector[loc, 0] += np.array(
                    [force_value*tck*L/2, 0, force_value*tck*L/2, 0])
                fc_type_dof = np.array([1, 0])
            elif force_dirc == 'fy':
                force_value_vector[loc, 0] += np.array(
                    [0, force_value*tck*L/2, 0, force_value*tck*L/2])
                fc_type_dof = np.array([0, 2])
    return force_value_vector, fc_type_dof


def force_surf(modelinfo, force_value, force_dirc, node_list_fc, dir_fc):
    if dir_fc == 'x':
        coord_fc = [3, 2]
    elif dir_fc == 'y':
        coord_fc = [1, 3]
    elif dir_fc == 'z':
        coord_fc = [1, 2]
    elmlist = np.array([0], dtype=int)
    for ii in range(len(node_list_fc)):
        elm2list = modelinfo['inci'][(np.asarray(
            np.where(modelinfo['inci'][:, 4:] == node_list_fc[ii])))[0][:], 0]
        elmlist = np.append(elmlist, elm2list)
    elmlist = np.unique(elmlist)
    elmlist = elmlist[1::][::]
    if modelinfo['elemid'][0] == 320:
        force_value_vector = np.zeros(
            (modelinfo['nodedof'][0]*len(node_list_fc), 1))
        for i in range(len(elmlist)):
            noi = int(modelinfo['inci'][int(elmlist[i]-1), 4])
            noj = int(modelinfo['inci'][int(elmlist[i]-1), 5])
            nok = int(modelinfo['inci'][int(elmlist[i]-1), 6])
            nol = int(modelinfo['inci'][int(elmlist[i]-1), 7])
            nom = int(modelinfo['inci'][int(elmlist[i]-1), 8])
            non = int(modelinfo['inci'][int(elmlist[i]-1), 9])
            noo = int(modelinfo['inci'][int(elmlist[i]-1), 10])
            nop = int(modelinfo['inci'][int(elmlist[i]-1), 11])
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
            coord_x_no1 = modelinfo['coord'][no1-1, 1]
            coord_y_no1 = modelinfo['coord'][no1-1, 2]
            coord_z_no1 = modelinfo['coord'][no1-1, 3]
            coord_x_no2 = modelinfo['coord'][no2-1, 1]
            coord_y_no2 = modelinfo['coord'][no2-1, 2]
            coord_z_no2 = modelinfo['coord'][no2-1, 3]
            coord_x_no3 = modelinfo['coord'][no3-1, 1]
            coord_y_no3 = modelinfo['coord'][no3-1, 2]
            coord_z_no3 = modelinfo['coord'][no3-1, 3]
            coord_x_no4 = modelinfo['coord'][no4-1, 1]
            coord_y_no4 = modelinfo['coord'][no4-1, 2]
            coord_z_no4 = modelinfo['coord'][no4-1, 3]
            poly = np.array([[coord_x_no1, coord_y_no1, coord_z_no1],
                             [coord_x_no2, coord_y_no2, coord_z_no2],
                             [coord_x_no3, coord_y_no3, coord_z_no3],
                             [coord_x_no4, coord_y_no4, coord_z_no4]])
            A = poly_area(poly)
            no1dof = np.where(no1 == node_list_fc)[0][0]
            no2dof = np.where(no2 == node_list_fc)[0][0]
            no3dof = np.where(no3 == node_list_fc)[0][0]
            no4dof = np.where(no4 == node_list_fc)[0][0]
            loc = np.array([modelinfo['nodedof'][0]*no1dof, modelinfo['nodedof'][0]*no1dof+1, modelinfo['nodedof'][0]*no1dof+2,
                            modelinfo['nodedof'][0]*no2dof, modelinfo['nodedof'][0] *
                            no2dof+1, modelinfo['nodedof'][0]*no2dof+2,
                            modelinfo['nodedof'][0]*no3dof, modelinfo['nodedof'][0] *
                            no3dof+1, modelinfo['nodedof'][0]*no3dof+2,
                            modelinfo['nodedof'][0]*no4dof, modelinfo['nodedof'][0]*no4dof+1, modelinfo['nodedof'][0]*no4dof+2, ])
            if force_dirc == 'fx':
                force_value_vector[loc, 0] += np.array(
                    [force_value*A/4, 0.0, 0.0, force_value*A/4, 0.0, 0.0, force_value*A/4, 0.0, 0.0, force_value*A/4, 0.0, 0.0])
                fc_type_dof = np.array([1, 0, 0])
            elif force_dirc == 'fy':
                force_value_vector[loc, 0] += np.array([0.0, force_value*A/4, 0.0, 0.0,
                                                       force_value*A/4, 0.0, 0.0, force_value*A/4, 0.0, 0.0, force_value*A/4, 0.0])
                fc_type_dof = np.array([0, 2, 0])
            elif force_dirc == 'fz':
                force_value_vector[loc, 0] += np.array([0.0, 0.0, force_value*A/4, 0.0, 0.0,
                                                       force_value*A/4, 0.0, 0.0, force_value*A/4, 0.0, 0.0, force_value*A/4])
                fc_type_dof = np.array([0, 0, 3])
    elif modelinfo['elemid'][0] == 310:
        force_value_vector = np.zeros(
            (modelinfo['nodedof'][0]*len(node_list_fc), 1))
        for i in range(len(elmlist)):
            noi = int(modelinfo['inci'][int(elmlist[i]-1), 4])
            noj = int(modelinfo['inci'][int(elmlist[i]-1), 5])
            nok = int(modelinfo['inci'][int(elmlist[i]-1), 6])
            nol = int(modelinfo['inci'][int(elmlist[i]-1), 7])
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
            # elif np.any([nok == node_list_fc[:]]):
            #     no1 = nok
            #     if np.any([nol == node_list_fc[:]]):
            #         no2 = nol
            else:
                continue
            coord_x_no1 = modelinfo['coord'][no1-1, 1]
            coord_y_no1 = modelinfo['coord'][no1-1, 2]
            coord_z_no1 = modelinfo['coord'][no1-1, 3]
            coord_x_no2 = modelinfo['coord'][no2-1, 1]
            coord_y_no2 = modelinfo['coord'][no2-1, 2]
            coord_z_no2 = modelinfo['coord'][no2-1, 3]
            coord_x_no3 = modelinfo['coord'][no3-1, 1]
            coord_y_no3 = modelinfo['coord'][no3-1, 2]
            coord_z_no3 = modelinfo['coord'][no3-1, 3]
            poly = np.array([[coord_x_no1, coord_y_no1, coord_z_no1],
                             [coord_x_no2, coord_y_no2, coord_z_no2],
                             [coord_x_no3, coord_y_no3, coord_z_no3]])
            A = poly_area(poly)
            no1dof = np.where(no1 == node_list_fc)[0][0]
            no2dof = np.where(no2 == node_list_fc)[0][0]
            no3dof = np.where(no3 == node_list_fc)[0][0]
            loc = np.array([modelinfo['nodedof'][0]*no1dof, modelinfo['nodedof'][0]*no1dof+1, modelinfo['nodedof'][0]*no1dof+2,
                            modelinfo['nodedof'][0]*no2dof, modelinfo['nodedof'][0] *
                            no2dof+1, modelinfo['nodedof'][0]*no2dof+2,
                            modelinfo['nodedof'][0]*no3dof, modelinfo['nodedof'][0]*no3dof+1, modelinfo['nodedof'][0]*no3dof+2])
            if force_dirc == 'fx':
                force_value_vector[loc, 0] += np.array(
                    [force_value*A/3, 0.0, 0.0, force_value*A/3, 0.0, 0.0, force_value*A/3, 0.0, 0.0])
                fc_type_dof = np.array([1, 0, 0])
            elif force_dirc == 'fy':
                force_value_vector[loc, 0] += np.array(
                    [0.0, force_value*A/3, 0.0, 0.0, force_value*A/3, 0.0, 0.0, force_value*A/3, 0.0])
                fc_type_dof = np.array([0, 2, 0])
            elif force_dirc == 'fz':
                force_value_vector[loc, 0] += np.array(
                    [0.0, 0.0, force_value*A/3, 0.0, 0.0, force_value*A/3, 0.0, 0.0, force_value*A/3])
                fc_type_dof = np.array([0, 0, 3])
    return force_value_vector, fc_type_dof


def force_beam(modelinfo, force_value, force_dirc, dir_fc, node_list_fc, line_fc):

    elmlist = [None]
    for ii in range(len(node_list_fc)):
        elm2list = modelinfo['inci'][(np.asarray(
            np.where(modelinfo['inci'][:, 4:] == node_list_fc[ii])))[0][:], 0]
        elmlist.extend(elm2list)
    elmlist = elmlist[1::][::]
    elmlist = np.unique(elmlist)
    if dir_fc == 'x':
        coord_fc = 1
    elif dir_fc == 'y':
        coord_fc = 2
    elif dir_fc == 'z':
        coord_fc = 3
    if modelinfo['elemid'][0] == 130:
        force_value_vector = np.zeros(
            (modelinfo['nodedof'][0]*len(node_list_fc), 1))
        for i in range(len(elmlist)):
            no1 = int(modelinfo['inci'][int(elmlist[i]-1), 4])
            no2 = int(modelinfo['inci'][int(elmlist[i]-1), 5])
            no1dof = np.asarray(np.where(no1 == node_list_fc))
            no2dof = np.asarray(np.where(no2 == node_list_fc))
            if no1dof.size != 0:
                if no2dof.size != 0:
                    L = np.sqrt((modelinfo['coord'][no1-1, 1] - modelinfo['coord'][no2-1, 1])**2 + (modelinfo['coord'][no1-1, 2] - modelinfo['coord'][no2-1, 2])**2 +
                                (modelinfo['coord'][no1-1, 3] - modelinfo['coord'][no2-1, 3])**2)

                    loc = np.array([modelinfo['nodedof'][0]*no1dof[0][0], modelinfo['nodedof'][0]*no1dof[0]
                                   [0]+1, modelinfo['nodedof'][0]*no2dof[0][0], modelinfo['nodedof'][0]*no2dof[0][0]+1])
                    force_value_vector[loc, 0] += np.array(
                        [force_value*L/2, (force_value*L**2)/12, force_value*L/2, (-1*force_value*L**2)/12])
                    fc_type_dof = np.array([2, 6])
    elif modelinfo['elemid'][0] == 140:
        force_value_vector = np.zeros(
            (modelinfo['nodedof'][0]*len(node_list_fc), 1))
        for i in range(len(elmlist)):
            no1 = int(modelinfo['inci'][int(elmlist[i]-1), 4])
            no2 = int(modelinfo['inci'][int(elmlist[i]-1), 5])
            no1dof = np.asarray(np.where(no1 == node_list_fc))
            no2dof = np.asarray(np.where(no2 == node_list_fc))
            if no1dof.size != 0:
                if no2dof.size != 0:
                    L = np.sqrt((modelinfo['coord'][no1-1, 1] - modelinfo['coord'][no2-1, 1])**2 + (modelinfo['coord'][no1-1, 2] - modelinfo['coord'][no2-1, 2])**2 +
                                (modelinfo['coord'][no1-1, 3] - modelinfo['coord'][no2-1, 3])**2)

                    loc = np.array([modelinfo['nodedof'][0]*no1dof[0][0], modelinfo['nodedof'][0]*no1dof[0][0]+1, modelinfo['nodedof'][0]*no1dof[0]
                                   [0]+2, modelinfo['nodedof'][0]*no2dof[0][0], modelinfo['nodedof'][0]*no2dof[0][0]+1, modelinfo['nodedof'][0]*no2dof[0][0]+2])
                    if force_dirc == 'fx':
                        force_value_vector[loc, 0] += np.array(
                            [(force_value*L**2)/2, 0.0, 0.0, (force_value*L**2)/2, 0.0, 0.0])
                        fc_type_dof = np.array([1, 0, 0])
                    elif force_dirc == 'fy':
                        force_value_vector[loc, 0] += np.array(
                            [0, force_value*L/2, (force_value*L**2)/12, 0, force_value*L/2, (-1*force_value*L**2)/12])
                        fc_type_dof = np.array([0, 2, 6])
    elif modelinfo['elemid'][0] == 141:
        force_value_vector = np.zeros(
            (modelinfo['nodedof'][0]*len(node_list_fc), 1))
        for i in range(len(elmlist)):
            no1 = int(modelinfo['inci'][int(elmlist[i]-1), 4])
            no2 = int(modelinfo['inci'][int(elmlist[i]-1), 5])
            no1dof = np.asarray(np.where(no1 == node_list_fc))
            no2dof = np.asarray(np.where(no2 == node_list_fc))
            if no1dof.size != 0:
                if no2dof.size != 0:
                    L = np.sqrt((modelinfo['coord'][no1-1, 1] - modelinfo['coord'][no2-1, 1])**2 + (modelinfo['coord'][no1-1, 2] - modelinfo['coord'][no2-1, 2])**2 +
                                (modelinfo['coord'][no1-1, 3] - modelinfo['coord'][no2-1, 3])**2)

                    loc = np.array([modelinfo['nodedof'][0]*no1dof[0][0], modelinfo['nodedof'][0]*no1dof[0][0]+1, modelinfo['nodedof'][0]*no1dof[0][0]+2, modelinfo['nodedof'][0]*no1dof[0][0]+3, modelinfo['nodedof'][0]*no1dof[0][0]+4, modelinfo['nodedof'][0]*no1dof[0][0]+5,
                                   modelinfo['nodedof'][0]*no2dof[0][0], modelinfo['nodedof'][0]*no2dof[0][0]+1, modelinfo['nodedof'][0]*no2dof[0][0]+2, modelinfo['nodedof'][0]*no2dof[0][0]+3, modelinfo['nodedof'][0]*no2dof[0][0]+4, modelinfo['nodedof'][0]*no2dof[0][0]+5])
                    if force_dirc == 'fx':
                        force_value_vector[loc, 0] += np.array(
                            [(force_value*L**2)/2, 0, 0, 0, 0, 0, (force_value*L**2)/2, 0, 0, 0, 0, 0])
                        fc_type_dof = np.array([1, 0, 0, 0, 0, 0])
                    elif force_dirc == 'fy':
                        force_value_vector[loc, 0] += np.array([0, force_value*L/2, 0, (
                            force_value*L**2)/12, 0, 0, 0, force_value*L/2, 0, (-1*force_value*L**2)/12, 0, 0])
                        fc_type_dof = np.array([0, 2, 0, 4, 0, 0])
                    elif force_dirc == 'fz':
                        force_value_vector[loc, 0] += np.array([0, 0, force_value*L/2, 0, (
                            force_value*L**2)/12, 0, 0, 0, force_value*L/2, 0, (-1*force_value*L**2)/120, 0])
                        fc_type_dof = np.array([0, 0, 3, 0, 5, 0])
    return force_value_vector, fc_type_dof


def poly_area(poly):
    if len(poly) < 3:  # not a plane - no area
        return 0
    total = [0, 0, 0]
    N = len(poly)
    for i in range(N):
        vi1 = poly[i]
        vi2 = poly[(i+1) % N]
        prod = np.cross(vi1, vi2)
        total[0] += prod[0]
        total[1] += prod[1]
        total[2] += prod[2]
    result = np.dot(total, unit_normal(poly[0], poly[1], poly[2]))
    return abs(result/2)


def unit_normal(a, b, c):
    x = np.linalg.det([[1, a[1], a[2]],
                       [1, b[1], b[2]],
                       [1, c[1], c[2]]])
    y = np.linalg.det([[a[0], 1, a[2]],
                       [b[0], 1, b[2]],
                       [c[0], 1, c[2]]])
    z = np.linalg.det([[a[0], a[1], 1],
                       [b[0], b[1], 1],
                       [c[0], c[1], 1]])
    magnitude = (x**2 + y**2 + z**2)**.5
    return (x/magnitude, y/magnitude, z/magnitude)
