#!/usr/bin/env python
__doc__ = """
get nodes
"""
import numpy as np


def search_edgex(edge_coordX: float, coord: np.ndarray, erro: float):
    """serch. node  on x dir. edge

    Arguments:
        edge_coordX:float   -- number coord in x dir.
        coord :np.array     -- nodes coordinates list in mesh
        erro:float          -- erro to conver.

    Returns:
        node                -- node loc.
    """
    dif = abs(edge_coordX * np.ones_like(coord[:, 1]) - coord[:, 1])
    node_posi = (np.asarray(np.where(dif < erro)))[0][:]
    node = coord[node_posi, 0].astype(int)
    return node


def search_edgey(edge_coordY: float, coord: np.ndarray, erro: float):
    """serch. node  on y dir. edge

    Arguments:
        edge_coordY:float   -- number coord in y dir.
        coord :np.array     -- nodes coordinates list in mesh
        erro:float          -- erro to conver.

    Returns:
        node                -- node loc.
    """
    dif = abs(edge_coordY * np.ones_like(coord[:, 2]) - coord[:, 2])
    node_posi = (np.asarray(np.where(dif < erro)))[0][:]
    node = coord[node_posi, 0].astype(int)
    return node


def search_edgez(edge_coordZ: float, coord: np.ndarray, erro: float):
    """serch. node  on z dir. edge

    Arguments:
        edge_coordY:float   -- number coord in z dir.
        coord :np.array     -- nodes coordinates list in mesh
        erro:float          -- erro to conver.

    Returns:
        node                -- node loc.
    """
    dif = abs(edge_coordZ * np.ones_like(coord[:, 3]) - coord[:, 3])
    node_posi = (np.asarray(np.where(dif < erro)))[0][:]
    node = coord[node_posi, 0].astype(int)
    return node


def search_surfxy(orthg_coordZ: float, coord: np.ndarray, erro: float):
    """serch. node on z dir. surf

    Arguments:
        orthg_coordZ:float  -- number coord in z dir.
        coord :np.array     -- nodes coordinates list in mesh
        erro:float          -- erro to conver.

    Returns:
        node                -- node loc.
    """
    dif = abs(orthg_coordZ * np.ones_like(coord[:, 3]) - coord[:, 3])
    node_posi = (np.asarray(np.where(dif < erro)))[0][:]
    node = coord[node_posi, 0].astype(int)
    return node


def search_surfyz(orthg_coordX: float, coord: np.ndarray, erro: float):
    """serch. node on x dir. surf

    Arguments:
        orthg_coordX:float  -- number coord in x dir.
        coord :np.array     -- nodes coordinates list in mesh
        erro:float          -- erro to conver.

    Returns:
        node                -- node loc.
    """
    dif = abs(orthg_coordX * np.ones_like(coord[:, 1]) - coord[:, 1])
    node_posi = (np.asarray(np.where(dif < erro)))[0][:]
    node = coord[node_posi, 0].astype(int)
    return node


def search_surfzx(orthg_coordY: float, coord: np.ndarray, erro: float):
    """serch. node ony dir. surf

    Arguments:
        orthg_coordY:float  -- number coord in y dir.
        coord :np.array     -- nodes coordinates list in mesh
        erro:float          -- erro to conver.

    Returns:
        node                -- node loc.
    """
    dif = abs(orthg_coordY * np.ones_like(coord[:, 2]) - coord[:, 2])
    node_posi = (np.asarray(np.where(dif < erro)))[0][:]
    node = coord[node_posi, 0].astype(int)
    return node


def search_nodexyz(node_coordX: float, node_coordY: float, node_coordZ: float, coord: np.ndarray, erro: float):
    """serch. node on coord mesh

    Arguments:
        node_coordX:float  -- number coord in x dir.
        node_coordY:float  -- number coord in y dir.
        node_coordZ:float  -- number coord in z dir.
        coord :np.array     -- nodes coordinates list in mesh
        erro:float          -- erro to conver.

    Returns:
        node                -- node loc.
    """
    difx = abs(node_coordX * np.ones_like(coord[:, 1]) - coord[:, 1])
    dify = abs(node_coordY * np.ones_like(coord[:, 2]) - coord[:, 2])
    difz = abs(node_coordZ * np.ones_like(coord[:, 3]) - coord[:, 3])
    node_posi = (np.asarray(np.where((difx < erro) & (dify < erro) & (difz < erro))))[
        0
    ][:]
    node = coord[node_posi, 0].astype(int)
    return node


def nodes_from_regions(regionlist: dict):
    """nodes from regions tag list

    Arguments:
        regionlist:dict    -- regions from tag list (gmsh mesh only)

    Returns:
        regions:dict
    """
    pointlist = regionlist["point"]
    points = [[None] * 2]
    for pp in range(len(pointlist)):
        points.append([pp + 1, np.unique(np.array(pointlist[str(pp + 1)]["nodes"]))])
    edgelist = regionlist["edge"]
    edges = [[None] * 2]
    for ee in range(len(edgelist)):
        edges.append([ee + 1, np.unique(np.array(edgelist[str(ee + 1)]["nodes"]))])
    surflist = regionlist["surf"]
    surfs = [[None] * 2]
    for ss in range(len(surflist)):
        surfs.append([ss + 1, np.unique(np.array(surflist[str(ss + 1)]["nodes"]))])
    points = points[1::][::]
    edges = edges[1::][::]
    surfs = surfs[1::][::]
    regions = [[None] * 2]
    regions.append(["P", points])
    regions.append(["E", edges])
    regions.append(["S", surfs])
    regions = regions[1::][::]
    return regions
