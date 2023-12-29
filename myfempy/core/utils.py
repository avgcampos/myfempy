import cython
import jax.numpy as jaxnp
import numpy as np
import scipy.sparse as sp

#==============================================================================
#                                      UTILS
#==============================================================================

def elem2nodes_conec(nnode: cython.int, nelem: cython.int, dofe: cython.int, inci):
    """
    Average Nodes Calculator version 2
    """
    # ith: cython.int[nelem * (dofe * dofe)]
    # jth: cython.int[nelem * (dofe * dofe)]
    # val: cython.double[nelem * (dofe * dofe)]

    ith = np.zeros((nelem * (dofe * dofe)), dtype=int)
    jth = np.zeros((nelem * (dofe * dofe)), dtype=int)
    val = np.zeros((nelem * (dofe * dofe)), dtype=int)
    q0 = 0
    for i in range(nnode):
        # elmlist = inci[(np.asarray(jaxnp.where(inci[:, 4:] == i + 1)))[0][:], 0]
        elmlist = inci[(np.asarray(np.where(inci[:, 4:] == i + 1)))[0][:], 0]
        q1 = elmlist.size
        ith[q0 : q1 + q0] = i
        jth[q0 : q1 + q0] = elmlist - 1
        val[q0 : q1 + q0] = elmlist
        q0 = q1 + q0
    S = sp.csc_matrix((val, (ith, jth)), shape=(nnode, nelem))
    return S


def results_average(results_elm, nnode: cython.int, nelem: cython.int, dofe: cython.int, inci):
    """_summary_

    Arguments:
        results_elm -- _description_

    Returns:
        _description_
    """
    # results_avr: cython.double[nnode]

    S = elem2nodes_conec(nnode, nelem, dofe, inci)
    results_avr = np.zeros((nnode), dtype=float)
    results_avr = [np.mean(results_elm[(S[mm, :].nonzero())[1]]) for mm in range(nnode)]

    # for mm in range(nnode):
    #     results_avr[mm] = np.mean(results_elm[(S[mm, :].nonzero())[1]])
    return results_avr


def poly_area(poly):
    if len(poly) < 3:  # not a plane - no area
        return 0
    total = [0, 0, 0]
    N = len(poly)
    for i in range(N):
        vi1 = poly[i]
        vi2 = poly[(i + 1) % N]
        prod = np.cross(vi1, vi2)
        total[0] += prod[0]
        total[1] += prod[1]
        total[2] += prod[2]
    result = np.dot(total, unit_normal(poly[0], poly[1], poly[2]))
    return abs(result / 2)


def unit_normal(a, b, c):
    x = np.linalg.det([[1, a[1], a[2]], [1, b[1], b[2]], [1, c[1], c[2]]])
    y = np.linalg.det([[a[0], 1, a[2]], [b[0], 1, b[2]], [c[0], 1, c[2]]])
    z = np.linalg.det([[a[0], a[1], 1], [b[0], b[1], 1], [c[0], c[1], 1]])
    magnitude = (x**2 + y**2 + z**2) ** 0.5
    return (x / magnitude, y / magnitude, z / magnitude)


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
    edgelist = regionlist["line"]
    edges = [[None] * 2]
    for ee in range(len(edgelist)):
        edges.append([ee + 1, np.unique(np.array(edgelist[str(ee + 1)]["nodes"]))])
    surflist = regionlist["plane"]
    surfs = [[None] * 2]
    for ss in range(len(surflist)):
        surfs.append([ss + 1, np.unique(np.array(surflist[str(ss + 1)]["nodes"]))])
    points = points[1::][::]
    edges = edges[1::][::]
    surfs = surfs[1::][::]
    regions = [[None] * 2]
    regions.append(["POINT_", points])
    regions.append(["LINE_", edges])
    regions.append(["PLANE_", surfs])
    regions = regions[1::][::]
    return regions




#integracao numerica
def gauss_points(type, npp):
    if type == 'quad4':
        return gauss_points_quad4(npp)
    
    elif type == 'tria3':
        return gauss_points_tria3(npp)
    
    elif type == 'hexa8':
        return gauss_points_hexa8(npp)
    
    elif type == 'tetr4':
        return gauss_points_tetr4(npp)
    
    else:
        pass


def gauss_points_quad4(npp):
    if npp == 1:
        pt = [[0, 0]]
        wt = [1.0]
        return pt, wt

    elif npp == 4:
        pt = [[-0.5773502691896258, -0.5773502691896258],
            [-0.5773502691896258, 0.5773502691896258],
            [0.5773502691896258, -0.5773502691896258],
            [0.5773502691896258, 0.5773502691896258]]
        wt = [1.0, 1.0, 1.0, 1.0]
        return pt, wt

    elif npp == 9:
        pt = [[-0.7745966692414834, -0.7745966692414834],
            [0.0000000000000000, -0.774596669241483],
            [0.7745966692414834, -0.7745966692414834],
            [-0.774596669241483, 0.0000000000000000],
            [0.0000000000000000, 0.0000000000000000],
            [0.7745966692414834, 0.0000000000000000],
            [-0.7745966692414834, 0.7745966692414834],
            [0.0000000000000000, 0.7745966692414834],
            [0.7745966692414834, 0.7745966692414834],
            ]
        wt = [0.555555555555556,
            0.888888888888889,
            0.555555555555556,
            0.888888888888889,
            0.888888888888889,
            0.888888888888889,
            0.555555555555556,
            0.888888888888889,
            0.555555555555556]
        return pt, wt


def gauss_points_tria3(npp):
    if npp == 1:
        pt = [[0.3333333333333333, 0.3333333333333333]]
        wt = [1.0]
        return pt, wt

    elif npp == 3:
        pt = [[0.1666666666666666, 0.1666666666666666],
              [0.6666666666666666, 0.1666666666666666],
              [0.1666666666666666, 0.6666666666666666]]
        wt = [0.3333333333333, 0.3333333333333, 0.3333333333333]
        return pt, wt

    elif npp == 7:
        pt = [[-0.7745966692414834, -0.7745966692414834],
            [0.0000000000000000, -0.774596669241483],
            [0.7745966692414834, -0.7745966692414834],
            [-0.774596669241483, 0.0000000000000000],
            [0.0000000000000000, 0.0000000000000000],
            [0.7745966692414834, 0.0000000000000000],
            [-0.7745966692414834, 0.7745966692414834],
            [0.0000000000000000, 0.7745966692414834],
            [0.7745966692414834, 0.7745966692414834],
            ]
        wt = [0.555555555555556,
            0.888888888888889,
            0.555555555555556,
            0.888888888888889,
            0.888888888888889,
            0.888888888888889,
            0.555555555555556,
            0.888888888888889,
            0.555555555555556]
        return pt, wt
    

def gauss_points_hexa8(npp):
    if npp == 1:
        pt = [[0, 0, 0]]
        wt = [1.0]
        return pt, wt

    elif npp == 8:
        pt = [[-0.5773502691896258, -0.5773502691896258, -0.5773502691896258],
              [0.5773502691896258, -0.5773502691896258, -0.5773502691896258],
              [0.5773502691896258, 0.5773502691896258, -0.5773502691896258],
              [-0.5773502691896258, 0.5773502691896258, -0.5773502691896258],
              [-0.5773502691896258, -0.5773502691896258, 0.5773502691896258],
              [0.5773502691896258, -0.5773502691896258, 0.5773502691896258],
              [0.5773502691896258, 0.5773502691896258, 0.5773502691896258],
              [-0.5773502691896258, 0.5773502691896258, 0.5773502691896258]]

        wt = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
        return pt, wt
    

def gauss_points_tetr4(npp):
    if npp == 1:
        pt = [[0.25, 0.25, 0.25]]
        wt = [1.0]
        return pt, wt

    elif npp == 4:
        pt = [[0.5854101966249685, 0.1381966011250105, 0.1381966011250105],
              [0.1381966011250105, 0.1381966011250105, 0.1381966011250105],
              [0.1381966011250105, 0.1381966011250105, 0.5854101966249685],
              [0.1381966011250105, 0.5854101966249685, 0.1381966011250105]]
        wt = [0.25, 0.25, 0.25, 0.25]
        return pt, wt

    elif npp == 5:
        pt = [[0.25, 0.25, 0.25],
              [0.5, 0.16666666666666666, 0.16666666666666666],
              [0.16666666666666666, 0.16666666666666666, 0.16666666666666666],
              [0.16666666666666666, 0.16666666666666666, 0.5],
              [0.16666666666666666, 0.5, 0.16666666666666666]]
        wt = [-0.8, 0.45, 0.45, 0.45, 0.45]
        return pt, wt