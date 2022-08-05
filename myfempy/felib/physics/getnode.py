#!/usr/bin/env python
__doc__ ="""
get nodes
"""
import numpy as np


def search_edgex(edge_coordX, coord, erro):
    dif = abs(edge_coordX*np.ones_like(coord[:, 1]) - coord[:, 1])
    node_posi = (np.asarray(np.where(dif < erro)))[0][:]
    node = coord[node_posi, 0].astype(int)
    return node


def search_edgey(edge_coordY, coord, erro):
    dif = abs(edge_coordY*np.ones_like(coord[:, 2]) - coord[:, 2])
    node_posi = (np.asarray(np.where(dif < erro)))[0][:]
    node = coord[node_posi, 0].astype(int)
    return node


def search_edgez(edge_coordZ, coord, erro):
    dif = abs(edge_coordZ*np.ones_like(coord[:, 3]) - coord[:, 3])
    node_posi = (np.asarray(np.where(dif < erro)))[0][:]
    node = coord[node_posi, 0].astype(int)
    return node


def search_surfxy(orthg_coordZ, coord, erro):
    dif = abs(orthg_coordZ*np.ones_like(coord[:, 3]) - coord[:, 3])
    node_posi = (np.asarray(np.where(dif < erro)))[0][:]
    node = coord[node_posi, 0].astype(int)
    return node


def search_surfyz(orthg_coordX, coord, erro):
    dif = abs(orthg_coordX*np.ones_like(coord[:, 1]) - coord[:, 1])
    node_posi = (np.asarray(np.where(dif < erro)))[0][:]
    node = coord[node_posi, 0].astype(int)
    return node


def search_surfzx(orthg_coordY, coord, erro):
    dif = abs(orthg_coordY*np.ones_like(coord[:, 2]) - coord[:, 2])
    node_posi = (np.asarray(np.where(dif < erro)))[0][:]
    node = coord[node_posi, 0].astype(int)
    return node


def search_nodexyz(node_coordX, node_coordY, node_coordZ, coord, erro):
    difx = abs(node_coordX*np.ones_like(coord[:, 1]) - coord[:, 1])
    dify = abs(node_coordY*np.ones_like(coord[:, 2]) - coord[:, 2])
    difz = abs(node_coordZ*np.ones_like(coord[:, 3]) - coord[:, 3])
    node_posi = (np.asarray(
        np.where((difx < erro) & (dify < erro) & (difz < erro))))[0][:]
    node = coord[node_posi, 0].astype(int)
    return node


def nodes_from_regions(regionlist):

    pointlist = regionlist['point']
    points = [[None]*2]
    for pp in range(len(pointlist)):
        points.append(
            [pp+1, np.unique(np.array(pointlist[str(pp+1)]['nodes']))])
    edgelist = regionlist['edge']
    edges = [[None]*2]
    for ee in range(len(edgelist)):
        edges.append([ee+1, np.unique(np.array(edgelist[str(ee+1)]['nodes']))])
    surflist = regionlist['surf']
    surfs = [[None]*2]
    for ss in range(len(surflist)):
        surfs.append([ss+1, np.unique(np.array(surflist[str(ss+1)]['nodes']))])
    points = points[1::][::]
    edges = edges[1::][::]
    surfs = surfs[1::][::]
    regions = [[None]*2]
    regions.append(['P', points])
    regions.append(['E', edges])
    regions.append(['S', surfs])
    regions = regions[1::][::]
    return regions
