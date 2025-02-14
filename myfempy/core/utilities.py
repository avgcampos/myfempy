
from numpy import (array, asarray, cross, dot, eye, float64, ix_, less, matmul,
                   mean, ones_like, sqrt, uint32, unique, where, zeros, empty)
from numpy.linalg import multi_dot
from scipy.linalg import block_diag, det, inv, kron
from scipy.sparse import csc_matrix

# from myfempy.core.util_cy import fast_dot


INT32 = uint32
FLT64 = float64

# ==============================================================================
#                               MYFEMPY UTTILITIES
# ==============================================================================

# def fastDOT(A, B):
#     R = empty((A.shape[0], B.shape[1]), dtype=float64)
#     fast_dot(A, B, R, 1)
#     return R

def addMatrix(A, A_add, idx):
    """
    addMatrix sum/add two matrix

    Arguments:
        A -- matrix numpy array
        A_add -- matrix numpy array to sum with A
        idx -- index dof relative

    Returns:
        matrix numpy array
    """

    A[ix_(idx, idx)] += A_add
    return A


def get2D_LocalVector(x):
    noix = x[0][0]
    noiy = x[1][0]
    nojx = x[2][0]
    nojy = x[3][0]
    L = sqrt((nojx - noix) ** 2 + (nojy - noiy) ** 2)
    s = (nojy - noiy) / L
    c = (nojx - noix) / L
    R = zeros((4, 4))
    R[0, 0] = c
    R[0, 1] = s
    R[1, 0] = -s
    R[1, 1] = c
    R[2, 2] = c
    R[2, 3] = s
    R[3, 2] = -s
    R[3, 3] = c
    xb = dot(R, x)
    return xb


def get3D_LocalVector(x, dim):
    noix = x[0][0]
    noiy = x[1][0]
    noiz = x[2][0]
    nojx = x[3][0]
    nojy = x[4][0]
    nojz = x[5][0]
    L = sqrt((nojx - noix) ** 2 + (nojy - noiy) ** 2 + (nojz - noiz) ** 2)
    if (noix == nojx) and (noiy == nojy):
        if nojz > noiz:
            lamb = array([[0, 0, 1.0], [0, 1.0, 0], [-1.0, 0, 0]])
        else:
            lamb = array([[0, 0, -1.0], [0, 1.0, 0], [1.0, 0, 0]])
    else:
        l = (nojx - noix) / L
        m = (nojy - noiy) / L
        n = (nojz - noiz) / L
        d = sqrt(l**2 + m**2)
        lamb = array([[l, m, n], [-m / d, l / d, 0.0], [-l * n / d, -m * n / d, d]])
    if dim == 1:
        R = block_diag(lamb)
        return dot(R, x)
    elif dim == 2:
        R = block_diag(lamb, lamb)
        return dot(R, x)
    elif dim == 3:
        R = block_diag(lamb, lamb, lamb)
        return dot(R, x)
    elif dim == 4:
        R = block_diag(lamb, lamb, lamb, lamb)
        return dot(R, x)
    elif dim == 5:
        R = block_diag(lamb, lamb, lamb, lamb, lamb)
        return dot(R, x)
    elif dim == 6:
        R = block_diag(lamb, lamb, lamb, lamb, lamb, lamb)
        return dot(R, x)
    else:
        return 0.0


def getRotational_Matrix(x, dim):
    noix = x[0][0]
    noiy = x[1][0]
    noiz = x[2][0]
    nojx = x[3][0]
    nojy = x[4][0]
    nojz = x[5][0]
    L = sqrt((nojx - noix) ** 2 + (nojy - noiy) ** 2 + (nojz - noiz) ** 2)
    if (noix == nojx) and (noiy == nojy):
        if nojz > noiz:
            lamb = array([[0, 0, 1.0], [0, 1.0, 0], [-1.0, 0, 0]])
        else:
            lamb = array([[0, 0, -1.0], [0, 1.0, 0], [1.0, 0, 0]])
    else:
        l = (nojx - noix) / L
        m = (nojy - noiy) / L
        n = (nojz - noiz) / L
        d = sqrt(l**2 + m**2)
        lamb = array([[l, m, n], [-m / d, l / d, 0.0], [-l * n / d, -m * n / d, d]])

    if dim == 1:
        return block_diag(lamb)
    elif dim == 2:
        return block_diag(lamb, lamb)
    elif dim == 3:
        return block_diag(lamb, lamb, lamb)
    elif dim == 4:
        return block_diag(lamb, lamb, lamb, lamb)
    elif dim == 5:
        return block_diag(lamb, lamb, lamb, lamb, lamb)
    elif dim == 6:
        return block_diag(lamb, lamb, lamb, lamb, lamb, lamb)
    else:
        return 0.0


def setSteps(steps):
    """
    setSteps steps setting

    Arguments:
        steps -- dict

    Returns:
        nsteps -- Number of steps [int]
    """
    start = steps["start"]
    end = steps["end"]
    substep = steps["step"]
    if (end - start) == 0:
        nsteps = int(end)
    else:
        nsteps = int((end - start) / substep)
    return nsteps


def elem2nodes_conec(nnode, nelem, dofe, inci):
    """
    Average Nodes Calculator version 2
    """
    ith = zeros((nelem * (dofe * dofe)), dtype=INT32)
    jth = zeros((nelem * (dofe * dofe)), dtype=INT32)
    val = zeros((nelem * (dofe * dofe)), dtype=INT32)
    q0 = 0
    for i in range(nnode):
        # elmlist = inci[(np.asarray(jaxnp.where(inci[:, 4:] == i + 1)))[0][:], 0]
        elmlist = inci[(asarray(where(inci[:, 4:] == i + 1)))[0][:], 0]
        q1 = elmlist.size
        ith[q0 : q1 + q0] = i
        jth[q0 : q1 + q0] = elmlist - 1
        val[q0 : q1 + q0] = elmlist
        q0 = q1 + q0
    S = csc_matrix((val, (ith, jth)), shape=(nnode, nelem), dtype=FLT64)
    return S


def results_average(results_elm, nnode, nelem, dofe, inci):
    """_summary_

    Arguments:
        results_elm -- _description_

    Returns:
        _description_
    """
    # results_avr: cython.double[nnode]

    S = elem2nodes_conec(nnode, nelem, dofe, inci)
    results_avr = zeros((nnode), dtype=FLT64)
    results_avr = [mean(results_elm[(S[mm, :].nonzero())[1]]) for mm in range(nnode)]

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
        prod = cross(vi1, vi2)
        total[0] += prod[0]
        total[1] += prod[1]
        total[2] += prod[2]
    result = dot(total, unit_normal(poly[0], poly[1], poly[2]))
    return abs(result / 2)


def unit_normal(a, b, c):
    x = det([[1, a[1], a[2]], [1, b[1], b[2]], [1, c[1], c[2]]])
    y = det([[a[0], 1, a[2]], [b[0], 1, b[2]], [c[0], 1, c[2]]])
    z = det([[a[0], a[1], 1], [b[0], b[1], 1], [c[0], c[1], 1]])
    magnitude = (x**2 + y**2 + z**2) ** 0.5
    return (x / magnitude, y / magnitude, z / magnitude)


def get_direction(coord, i, j):
    direc = zeros(3)
    rank = coord.shape[1]
    direc[:rank] = coord[i, :] - coord[j, :]
    return direc


def get_nodes_from_list(nodelist, coord, regions):
    tol = 1e-10
    nodes = []
    dir_fc = []
    # ----- SEEKERS WITH LOC -----
    if nodelist[0] == "lengthx":
        coord_0 = float(nodelist[2])
        coord_1 = float(nodelist[3])
        nodes = coord[
            where((coord[:, 1] >= coord_0) & (coord[:, 1] <= coord_1)),
            0,
        ][0]
        # nodesapply.append(nodes)
        dir_fc = "x"
        return nodes, dir_fc

    elif nodelist[0] == "lengthy":
        coord_0 = float(nodelist[2])
        coord_1 = float(nodelist[3])
        nodes = coord[
            where((coord[:, 2] >= coord_0) & (coord[:, 2] <= coord_1)),
            0,
        ][0]
        # nodesapply.append(nodes)
        dir_fc = "y"
        return nodes, dir_fc

    elif nodelist[0] == "lengthz":
        coord_0 = float(nodelist[2])
        coord_1 = float(nodelist[3])
        nodes = coord[
            where((coord[:, 3] >= coord_0) & (coord[:, 3] <= coord_1)),
            0,
        ][0]
        # nodesapply.append(nodes)
        dir_fc = "z"
        return nodes, dir_fc

    elif nodelist[0] == "edgex":
        edge_coordX = float(
            nodelist[1]
        )  # nodelist[1].astype(float)  # float(nodelist[1])
        if int(nodelist[2]) == 999:
            dir_fc = "x_y"
            coord_fc = (
                coord[
                    where(coord[:, 3] == float(nodelist[3])),
                    :,
                ]
            )[0]
        elif int(nodelist[3]) == 999:
            dir_fc = "x_z"
            coord_fc = (
                coord[
                    where(coord[:, 2] == float(nodelist[2])),
                    :,
                ]
            )[0]
        nodes = search_edgex(edge_coordX, coord_fc, tol)
        # nodesapply.append(nodes)
        return nodes, dir_fc

    elif nodelist[0] == "edgey":
        edge_coordY = float(nodelist[2])
        if float(nodelist[1]) == 999:
            dir_fc = "y_x"
            coord_fc = (
                coord[
                    where(coord[:, 3] == float(nodelist[3])),
                    :,
                ]
            )[0]
        elif float(nodelist[3]) == 999:
            dir_fc = "y_z"
            coord_fc = (
                coord[
                    where(coord[:, 1] == float(nodelist[1])),
                    :,
                ]
            )[0]
        nodes = search_edgey(edge_coordY, coord_fc, tol)
        # nodesapply.append(nodes)
        return nodes, dir_fc

    elif nodelist[0] == "edgez":
        edge_coordZ = float(nodelist[3])
        if float(nodelist[1]) == 999:
            dir_fc = "z_x"
            coord_fc = (
                coord[
                    where(coord[:, 2] == float(nodelist[2])),
                    :,
                ]
            )[0]

        elif float(nodelist[2]) == 999:
            dir_fc = "z_x"
            coord_fc = (
                coord[
                    where(coord[:, 1] == float(nodelist[1])),
                    :,
                ]
            )[0]
        nodes = search_edgez(edge_coordZ, coord_fc, tol)
        return nodes, dir_fc
        # nodesapply.append(nodes)

    elif nodelist[0] == "surfxy":
        orthg_coordZ = float(nodelist[3])
        nodes = search_surfxy(orthg_coordZ, coord, tol)
        # nodesapply.append(nodes)
        dir_fc = "z"
        return nodes, dir_fc

    elif nodelist[0] == "surfyz":
        orthg_coordX = float(nodelist[1])
        nodes = search_surfyz(orthg_coordX, coord, tol)
        # nodesapply.append(nodes)
        dir_fc = "x"
        return nodes, dir_fc

    elif nodelist[0] == "surfzx":
        orthg_coordY = float(nodelist[2])
        nodes = search_surfzx(orthg_coordY, coord, tol)
        # nodesapply.append(nodes)
        dir_fc = "y"
        return nodes, dir_fc

    elif nodelist[0] == "node":
        node_coordX = float(nodelist[1])
        node_coordY = float(nodelist[2])
        node_coordZ = float(nodelist[3])
        nodes = search_nodexyz(node_coordX, node_coordY, node_coordZ, coord, tol)
        # nodesapply.append(nodes)
        dir_fc = "x"
        return nodes, dir_fc

    # ----- SEEKERS WITH TAG -----
    elif nodelist[0] == "point":
        nodes = regions[0][1][int(nodelist[4]) - 1][1][:]
        dir_fc = "x"
        # nodesapply.append(nodes)
        return nodes, dir_fc

    elif nodelist[0] == "line":
        nodes = regions[1][1][int(nodelist[4]) - 1][1][:]
        dir_fc = "x"
        # nodesapply.append(nodes)
        return nodes, dir_fc

    elif nodelist[0] == "plane":
        nodes = regions[2][1][int(nodelist[4]) - 1][1][:]
        dir_fc = "x"
        # nodesapply.append(nodes)
        return nodes, dir_fc

    else:
        print("input erro: force_opt don't defined")
        return nodes, dir_fc


def get_elemen_from_nodelist(inci, node_list):
    elmlist = [None]
    for ii in range(len(node_list)):
        elm2list = inci[(asarray(where(inci[:, 4:] == node_list[ii])))[0][:], 0]
        elmlist.extend(elm2list)
    elmlist = elmlist[1::][::]
    elmlist = unique(elmlist)
    return elmlist


def search_edgex(edge_coordX, coord, erro):
    """serch. node  on x dir. edge

    Arguments:
        edge_coordX:float   -- number coord in x dir.
        coord :np.array     -- nodes coordinates list in mesh
        erro:float          -- erro to conver.

    Returns:
        node                -- node loc.
    """
    dif = abs(edge_coordX * ones_like(coord[:, 1]) - coord[:, 1])
    erroa = erro * ones_like(dif)
    node_posi = less(
        dif, erroa
    )  # array(where(dif < erroa)[0][:]).astype(INT32) #(asarray(where(dif < erroa)))[0][:]
    node = coord[node_posi, 0].astype(INT32)
    return node


def search_edgey(edge_coordY, coord, erro):
    """serch. node  on y dir. edge

    Arguments:
        edge_coordY:float   -- number coord in y dir.
        coord :np.array     -- nodes coordinates list in mesh
        erro:float          -- erro to conver.

    Returns:
        node                -- node loc.
    """
    dif = abs(edge_coordY * ones_like(coord[:, 2]) - coord[:, 2])
    erroa = erro * ones_like(dif)
    node_posi = less(dif, erroa)  # (asarray(where(dif < erroa)))[0][:]
    node = coord[node_posi, 0].astype(INT32)
    return node


def search_edgez(edge_coordZ, coord, erro):
    """serch. node  on z dir. edge

    Arguments:
        edge_coordY:float   -- number coord in z dir.
        coord :np.array     -- nodes coordinates list in mesh
        erro:float          -- erro to conver.

    Returns:
        node                -- node loc.
    """
    dif = abs(edge_coordZ * ones_like(coord[:, 3]) - coord[:, 3])
    erroa = erro * ones_like(dif)
    node_posi = less(dif, erroa)  # (asarray(where(dif < erroa)))[0][:]
    node = coord[node_posi, 0].astype(INT32)
    return node


def search_surfxy(orthg_coordZ, coord, erro):
    """serch. node on z dir. surf

    Arguments:
        orthg_coordZ:float  -- number coord in z dir.
        coord :np.array     -- nodes coordinates list in mesh
        erro:float          -- erro to conver.

    Returns:
        node                -- node loc.
    """
    dif = abs(orthg_coordZ * ones_like(coord[:, 3]) - coord[:, 3])
    erroa = erro * ones_like(dif)
    node_posi = less(dif, erroa)  # (asarray(where(dif < erroa)))[0][:]
    node = coord[node_posi, 0].astype(INT32)
    return node


def search_surfyz(orthg_coordX, coord, erro):
    """serch. node on x dir. surf

    Arguments:
        orthg_coordX:float  -- number coord in x dir.
        coord :np.array     -- nodes coordinates list in mesh
        erro:float          -- erro to conver.

    Returns:
        node                -- node loc.
    """
    dif = abs(orthg_coordX * ones_like(coord[:, 1]) - coord[:, 1])
    erroa = erro * ones_like(dif)
    node_posi = less(dif, erroa)  # (asarray(where(dif < erroa)))[0][:]
    node = coord[node_posi, 0].astype(INT32)
    return node


def search_surfzx(orthg_coordY, coord, erro):
    """serch. node ony dir. surf

    Arguments:
        orthg_coordY:float  -- number coord in y dir.
        coord :np.array     -- nodes coordinates list in mesh
        erro:float          -- erro to conver.

    Returns:
        node                -- node loc.
    """
    dif = abs(orthg_coordY * ones_like(coord[:, 2]) - coord[:, 2])
    erroa = erro * ones_like(dif)
    node_posi = less(dif, erroa)  # (asarray(where(dif < erroa)))[0][:]
    node = coord[node_posi, 0].astype(INT32)
    return node


def search_nodexyz(node_coordX, node_coordY, node_coordZ, coord, erro):
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
    difx = abs(node_coordX * ones_like(coord[:, 1]) - coord[:, 1])
    dify = abs(node_coordY * ones_like(coord[:, 2]) - coord[:, 2])
    difz = abs(node_coordZ * ones_like(coord[:, 3]) - coord[:, 3])
    erroax = erro * ones_like(difx)
    erroay = erro * ones_like(dify)
    erroaz = erro * ones_like(difz)
    node_posi = (asarray(where((difx < erroax) & (dify < erroay) & (difz < erroaz))))[
        0
    ][:]
    node = coord[node_posi, 0].astype(INT32)
    return node


def nodes_from_regions(regionlist):
    """nodes from regions tag list

    Arguments:
        regionlist:dict    -- regions from tag list (gmsh mesh only)

    Returns:
        regions:dict
    """
    pointlist = regionlist["point"]
    points = [[None] * 2]
    for pp in range(len(pointlist)):
        points.append([pp + 1, unique(array(pointlist[str(pp + 1)]["nodes"]))])
    edgelist = regionlist["line"]
    edges = [[None] * 2]
    for ee in range(len(edgelist)):
        edges.append([ee + 1, unique(array(edgelist[str(ee + 1)]["nodes"]))])
    surflist = regionlist["plane"]
    surfs = [[None] * 2]
    for ss in range(len(surflist)):
        surfs.append([ss + 1, unique(array(surflist[str(ss + 1)]["nodes"]))])
    points = points[1::][::]
    edges = edges[1::][::]
    surfs = surfs[1::][::]
    regions = [[None] * 2]
    regions.append(["POINT_", points])
    regions.append(["LINE_", edges])
    regions.append(["PLANE_", surfs])
    regions = regions[1::][::]
    return regions


# integracao numerica
def gauss_points(type, npp):

    if type == "line2":
        return __gauss_points_line(npp)

    elif type == "line3":
        return __gauss_points_line(npp)

    elif type == "tria3":
        return __gauss_points_tria(npp)

    elif type == "tria6":
        return __gauss_points_tria(npp)

    elif type == "quad4":
        return __gauss_points_quad(npp)

    elif type == "quad8":
        return __gauss_points_quad(npp)

    elif type == "hexa8":
        return __gauss_points_quad(npp)

    elif type == "tetr4":
        return __gauss_points_tetr(npp)

    else:
        pass


def __gauss_points_line(npp):
    if npp == 1:
        pt = array([0.000000000000000])
        wt = array([2.000000000000000])
        return pt, wt

    elif npp == 2:
        pt = array([-0.5773502691896258, 0.5773502691896258])
        wt = array([1.000000000000000, 1.000000000000000])
        return pt, wt

    elif npp == 4:
        pt = array([-0.86113631, -0.33998104, 0.33998104, 0.86113631])
        wt = array([0.34785485, 0.65214515, 0.65214515, 0.34785485])
        return pt, wt

    elif npp == 8:
        pt = array(
            [
                -0.96028986,
                -0.79666648,
                -0.52553241,
                -0.18343464,
                0.18343464,
                0.52553241,
                0.79666648,
                0.96028986,
            ]
        )
        wt = array(
            [
                0.10122854,
                0.22238103,
                0.31370665,
                0.36268378,
                0.36268378,
                0.31370665,
                0.22238103,
                0.10122854,
            ]
        )

        return pt, wt


def __gauss_points_quad(npp):
    if npp == 1:
        pt = array([0.000000000000000])
        wt = array([2.000000000000000])
        return pt, wt

    elif npp == 2:
        pt = array([-0.5773502691896258, 0.5773502691896258])
        wt = array([1.000000000000000, 1.000000000000000])
        return pt, wt

    elif npp == 3:
        pt = array([-0.774596669241483, 0.000000000000000, 0.774596669241483])
        wt = array([0.555555555555556, 0.888888888888889, 0.555555555555556])
        return pt, wt

    elif npp == 4:
        pt = array([-0.86113631, -0.33998104, 0.33998104, 0.86113631])
        wt = array([0.34785485, 0.65214515, 0.65214515, 0.34785485])
        return pt, wt

    elif npp == 8:
        pt = array(
            [
                -0.96028986,
                -0.79666648,
                -0.52553241,
                -0.18343464,
                0.18343464,
                0.52553241,
                0.79666648,
                0.96028986,
            ]
        )
        wt = array(
            [
                0.10122854,
                0.22238103,
                0.31370665,
                0.36268378,
                0.36268378,
                0.31370665,
                0.22238103,
                0.10122854,
            ]
        )
        return pt, wt

    elif npp == 9:
        pt = array(
            [
                -0.96816024,
                -0.83603111,
                -0.61337143,
                -0.32425342,
                0.0,
                0.32425342,
                0.61337143,
                0.83603111,
                0.96816024,
            ]
        )
        wt = array(
            [
                0.08127439,
                0.18064816,
                0.2606107,
                0.31234708,
                0.33023936,
                0.31234708,
                0.2606107,
                0.18064816,
                0.08127439,
            ]
        )
        return pt, wt


# def __gauss_points_hexa(npp):
#     if npp == 1:
#         pt = array([[0.000000000000000]])
#         wt = array([2.000000000000000])
#         return pt, wt

#     elif npp == 8:
#         pt = array(
#             [
#                 [-0.5773502691896258, -0.5773502691896258, -0.5773502691896258],
#                 [0.5773502691896258, -0.5773502691896258, -0.5773502691896258],
#                 [0.5773502691896258, 0.5773502691896258, -0.5773502691896258],
#                 [-0.5773502691896258, 0.5773502691896258, -0.5773502691896258],
#                 [-0.5773502691896258, -0.5773502691896258, 0.5773502691896258],
#                 [0.5773502691896258, -0.5773502691896258, 0.5773502691896258],
#                 [0.5773502691896258, 0.5773502691896258, 0.5773502691896258],
#                 [-0.5773502691896258, 0.5773502691896258, 0.5773502691896258],
#             ]
#         )

#         wt = array([1.000000000000000,
#                     1.000000000000000,
#                     1.000000000000000,
#                     1.000000000000000,
#                     1.000000000000000,
#                     1.0000000000000000,
#                     1.000000000000000,
#                     1.000000000000000])
#         return pt, wt


def __gauss_points_tria(npp):
    if npp == 1:
        pt = array([0.3333333333333333])
        wt = array([1.000000000000000])
        return pt, wt

    elif npp == 3:
        pt = array([0.6666666666666666, 0.1666666666666666, 0.1666666666666666])
        wt = array([0.3333333333333, 0.3333333333333, 0.3333333333333])
        return pt, wt

    elif npp == 7:
        pt = array(
            [
                0.1012865073235,
                0.79742698553531,
                0.1012865073235,
                0.4701420641051,
                0.4701420641051,
                0.597158717898,
                0.3333333333333,
            ]
        )
        wt = array(
            [
                0.1259391805448,
                0.1259391805448,
                0.1259391805448,
                0.1323941527885,
                0.1323941527885,
                0.1323941527885,
                0.2250000000000,
            ]
        )
        return pt, wt


def __gauss_points_tetr(npp):
    if npp == 1:
        pt = array([0.3333333333333333])
        wt = array([0.166666666666667])
        return pt, wt
    elif npp == 4:
        pt = array(
            [
                0.5854101966249685,
                0.1381966011250105,
                0.1381966011250105,
                0.1381966011250105,
            ]
        )
        wt = array(
            [
                0.25000000000000000,
                0.25000000000000000,
                0.25000000000000000,
                0.25000000000000000,
            ]
        )
        return pt, wt
    elif npp == 5:
        pt = array(
            [
                0.25000000000000000,
                0.50000000000000000,
                0.166666666666667,
                0.166666666666667,
                0.166666666666667,
            ]
        )
        wt = array(
            [
                -0.8000000000000000,
                0.45000000000000000,
                0.45000000000000000,
                0.45000000000000000,
                0.45000000000000000,
            ]
        )
        return pt, wt


# antigo

# def __gauss_points_quad(npp):
#     if npp == 1:
#         pt = array([[0.000000000000000, 0.000000000000000]])
#         wt = array([2.000000000000000])
#         return pt, wt

#     elif npp == 2:
#         pt = array(
#             [
#                 [-0.5773502691896258, -0.5773502691896258],
#                 [0.5773502691896258, 0.5773502691896258],
#             ]
#         )
#         wt = array([1.000000000000000,
#                     1.000000000000000])
#         return pt, wt

#     elif npp == 4:
#         pt = array(
#             [
#                 [-0.5773502691896258, -0.5773502691896258],
#                 [-0.5773502691896258, 0.5773502691896258],
#                 [0.5773502691896258, -0.5773502691896258],
#                 [0.5773502691896258, 0.5773502691896258],
#             ]
#         )
#         wt = array([1.000000000000000,
#                     1.000000000000000,
#                     1.000000000000000,
#                     1.000000000000000])
#         return pt, wt

#     elif npp == 9:
#         pt = array(
#             [
#                 [-0.774596669241483, -0.774596669241483],
#                 [0.000000000000000, -0.77459666924148],
#                 [0.774596669241483, -0.774596669241483],
#                 [-0.774596669241483, 0.000000000000000],
#                 [0.000000000000000, 0.000000000000000],
#                 [0.774596669241483, 0.000000000000000],
#                 [-0.774596669241483, 0.774596669241483],
#                 [0.000000000000000, 0.774596669241483],
#                 [0.774596669241483, 0.774596669241483],
#             ]
#         )
#         wt = array(
#             [
#                 0.555555555555556,
#                 0.888888888888889,
#                 0.555555555555556,
#                 0.888888888888889,
#                 0.888888888888889,
#                 0.888888888888889,
#                 0.555555555555556,
#                 0.888888888888889,
#                 0.555555555555556,
#             ]
#         )
#         return pt, wt


# def __gauss_points_tria(npp):
#     if npp == 1:
#         pt = array([[0.3333333333333333, 0.3333333333333333]])
#         wt = array([1.000000000000000])
#         return pt, wt

#     elif npp == 2:
#         pt = array([[0.5000000000000000, 0.5000000000000000], [0.5000000000000000, 0.5000000000000000]])
#         wt = array([1.000000000000000,
#                     1.000000000000000])
#         return pt, wt

#     elif npp == 3:
#         pt = array(
#             [
#                 [0.1666666666666666, 0.1666666666666666],
#                 [0.6666666666666666, 0.1666666666666666],
#                 [0.1666666666666666, 0.6666666666666666],
#             ]
#         )
#         wt = array([0.3333333333333,
#                     0.3333333333333,
#                     0.3333333333333])
#         return pt, wt

#     elif npp == 7:
#         pt = array(
#             [
#                 [-0.7745966692414834, -0.7745966692414834],
#                 [0.0000000000000000, -0.774596669241483],
#                 [0.7745966692414834, -0.7745966692414834],
#                 [-0.774596669241483, 0.0000000000000000],
#                 [0.0000000000000000, 0.0000000000000000],
#                 [0.7745966692414834, 0.0000000000000000],
#                 [-0.7745966692414834, 0.7745966692414834],
#                 [0.0000000000000000, 0.7745966692414834],
#                 [0.7745966692414834, 0.7745966692414834],
#             ]
#         )
#         wt = array(
#             [
#                 0.555555555555556,
#                 0.888888888888889,
#                 0.555555555555556,
#                 0.888888888888889,
#                 0.888888888888889,
#                 0.888888888888889,
#                 0.555555555555556,
#                 0.888888888888889,
#                 0.555555555555556,
#             ]
#         )
#         return pt, wt


#     elif npp == 5:
#         pt = array(
#             [
#                 [0.25000000000000000, 0.25000000000000000, 0.25000000000000000],
#                 [0.5000000000000000, 0.16666666666666666, 0.16666666666666666],
#                 [0.16666666666666666, 0.16666666666666666, 0.16666666666666666],
#                 [0.16666666666666666, 0.16666666666666666, 0.5000000000000000],
#                 [0.16666666666666666, 0.5000000000000000, 0.16666666666666666],
#             ]
#         )
#         wt = array([-0.8000000000000000,
#                     0.45000000000000000,
#                     0.45000000000000000,
#                     0.45000000000000000,
#                     0.45000000000000000])
#         return pt, wt

