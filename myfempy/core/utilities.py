# from os import environ
# environ['OMP_NUM_THREADS'] = '3'

from numpy import array, zeros, eye, dot, matmul, asarray, where, mean, cross, ones_like, unique, less, ix_, uint32, float64
from numpy.linalg import multi_dot
from scipy.linalg import inv, kron, det
from scipy.sparse import csc_matrix
INT32 = uint32
FLT64 = float64

#==============================================================================
#                               MYFEMPY UTTILITIES
#==============================================================================
# def inverse(A):
#     invA = inv(A, overwrite_a=True, check_finite=False)
#     return invA

# def kronProd(A, B):
#     kronAB = kron(A, B)
#     return kronAB

# def determinant(A):
#     detA = det(A, overwrite_a=True, check_finite=False)
#     return detA

# def dotProd(A, B):
#     dotAB = matmul(A, B) #dot(A, B)tAB
#     return dotAB

# def dotdotProd(A, B, C):
#     mdotABC = multi_dot([A, B, C])
#     return mdotABC

# def getZerosArray(m, n, type):
#     zeros_array = zeros((m, n), dtype=type)
#     return zeros_array

# def getNewArray(array_list, type):
#     new_array = array(array_list, dtype=type)
#     return new_array

# def getEyeMatrix(n, type):
#     eyeM = eye(n, dtype=type)
#     return eyeM

def determinant_dim2(A):
    detA = A[0]*A[3]-A[1]*A[2]
    return detA

def inverse_dim2(A):
    invA = 1/(A[0]*A[3]-A[1]*A[2])*array([[A[3], -A[1]], [-A[2], A[0]]])
    return invA


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
    # ith: cython.int[nelem * (dofe * dofe)]
    # jth: cython.int[nelem * (dofe * dofe)]
    # val: cython.double[nelem * (dofe * dofe)]

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



def get_nodes_from_list(nodelist, coord, regions):
    tol = 1e-10
    nodes=[]
    dir_fc=[]
    # ----- SEEKERS WITH LOC -----
    if nodelist[0] == "lengthx":
        coord_0 = float(nodelist[2])
        coord_1 = float(nodelist[3])
        nodes = coord[where((coord[:, 1] >= coord_0)&(coord[:, 1] <= coord_1)), 0,][0]
        # nodesapply.append(nodes)
        dir_fc = "x"
        return nodes, dir_fc
    
    elif nodelist[0] == "lengthy":
        coord_0 = float(nodelist[2])
        coord_1 = float(nodelist[3])
        nodes = coord[where((coord[:, 2] >= coord_0)&(coord[:, 2] <= coord_1)), 0,][0]
        # nodesapply.append(nodes)
        dir_fc = "y"
        return nodes, dir_fc
    
    elif nodelist[0] == "lengthz":
        coord_0 = float(nodelist[2])
        coord_1 = float(nodelist[3])
        nodes = coord[where((coord[:, 3] >= coord_0)&(coord[:, 3] <= coord_1)), 0,][0]
        # nodesapply.append(nodes)
        dir_fc = "z"
        return nodes, dir_fc
    
    elif nodelist[0] == "edgex":
        edge_coordX = nodelist[1].astype(float) #float(nodelist[1])
        if int(nodelist[2]) == 999:
            dir_fc = "x_y"
            coord_fc = (coord[where(coord[:, 3] == float(nodelist[3])),:,])[0]
        elif int(nodelist[3]) == 999:
            dir_fc = "x_z"
            coord_fc = (coord[where(coord[:, 2] == float(nodelist[2])),:,])[0]
        nodes = search_edgex(edge_coordX, coord_fc, tol)
        # nodesapply.append(nodes)
        return nodes, dir_fc
    
    elif nodelist[0] == "edgey":
        edge_coordY = float(nodelist[2])
        if float(nodelist[1]) == 999:
            dir_fc = "y_x"
            coord_fc = (coord[where(coord[:, 3] == float(nodelist[3])),:,])[0]
        elif float(nodelist[3]) == 999:
            dir_fc = "y_z"
            coord_fc = (coord[where(coord[:, 1] == float(nodelist[1])),:,])[0]
        nodes = search_edgey(edge_coordY, coord_fc, tol)
        # nodesapply.append(nodes)
        return nodes, dir_fc
    
    elif nodelist[0] == "edgez":
        edge_coordZ = float(nodelist[3])
        if float(nodelist[1]) == 999:
            dir_fc = "z_x"
            coord_fc = (coord[where(coord[:, 2] == float(nodelist[2])),:,])[0]
            
        elif float(nodelist[2]) == 999:
            dir_fc = "z_x"
            coord_fc = (coord[where(coord[:, 1] == float(nodelist[1])), :,])[0]
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
        nodes = regions[0][1][int(nodelist[-1]) - 1][1][:]
        dir_fc = "x"
        # nodesapply.append(nodes)
        return nodes, dir_fc
    
    elif nodelist[0] == "line":
        nodes = regions[1][1][int(nodelist[-1]) - 1][1][:]
        dir_fc = "x"
        # nodesapply.append(nodes)
        return nodes, dir_fc
    
    elif nodelist[0] == "plane":
        nodes = regions[2][1][int(nodelist[-1]) - 1][1][:]
        dir_fc = "x"
        # nodesapply.append(nodes)
        return nodes, dir_fc
    
    else:
        print("input erro: force_opt don't defined")
        return nodes, dir_fc


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
    node_posi = less(dif, erroa) #array(where(dif < erroa)[0][:]).astype(INT32) #(asarray(where(dif < erroa)))[0][:]
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
    node_posi = less(dif, erroa) #(asarray(where(dif < erroa)))[0][:]
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
    node_posi = less(dif, erroa) #(asarray(where(dif < erroa)))[0][:]
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
    node_posi = less(dif, erroa) #(asarray(where(dif < erroa)))[0][:]
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
    node_posi = less(dif, erroa) #(asarray(where(dif < erroa)))[0][:]
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
    node_posi = less(dif, erroa) #(asarray(where(dif < erroa)))[0][:]
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
    node_posi = (asarray(where((difx < erroax) & (dify < erroay) & (difz < erroaz))))[0][:]
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


#integracao numerica
def gauss_points(type, npp):
    if type == 'quad4':
        return __gauss_points_quad4(npp)
    
    elif type == 'tria3':
        return __gauss_points_tria3(npp)
    
    elif type == 'hexa8':
        return __gauss_points_hexa8(npp)
    
    elif type == 'tetr4':
        return __gauss_points_tetr4(npp)
    
    else:
        pass


def __gauss_points_quad4(npp):
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


def __gauss_points_tria3(npp):
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
    

def __gauss_points_tetr4(npp):
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
    

def __gauss_points_hexa8(npp):
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
    