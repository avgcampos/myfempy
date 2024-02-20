from os import environ
environ['OMP_NUM_THREADS'] = '3'

from numpy import array, zeros, eye, dot, asarray, where, mean, cross, ones_like, unique, less, uint32, float64
from scipy.linalg import inv, kron, det
from scipy.sparse import csc_matrix
cimport cython
INT32 = uint32
FLT64 = float64

#==============================================================================
#                               MYFEMPY UTTILITIES
#==============================================================================
@cython.exceptval(check=False)
@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
def inverse(double[:, :] A):
    cdef double [:, :] invA = inv(A)
    return asarray(invA)

@cython.exceptval(check=False)
@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
def kronProd(double[:, :] A, double[:, :] B):
    cdef double [:, :] kronAB = kron(A, B)
    return asarray(kronAB)

@cython.exceptval(check=False)
@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
def determinant(double[:, :] A):
    cdef double detA = det(A)
    return asarray(detA)

@cython.exceptval(check=False)
@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
def dotProd(double[:, :] A, double[:, :] B):
    cdef double [:, :] dotAB = dot(A, B)
    return asarray(dotAB)

@cython.exceptval(check=False)
@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
def getZerosArray(int m, int n, type):
    cdef double [:, :] A = zeros((m, n), dtype=type)
    return asarray(A)

@cython.exceptval(check=False)
@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
def getNewArray(list array_list, type):
    cdef double [:, :] A = array(array_list, dtype=type)
    return asarray(A)

@cython.exceptval(check=False)
@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
def getEyeMatrix(int n, type):
    cdef unsigned int [:, :] A = eye(n, dtype=type)
    return asarray(A)
    
@cython.exceptval(check=False)
@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
def elem2nodes_conec(int nnode, int nelem, int dofe, double[:, :] inci):
    cdef int i, node, q0, q1, nelemxdofexdofe
    cdef double[:, :] S
    nelemxdofexdofe = nelem * dofe * dofe
    cdef int[:, ::1] ith = zeros((nelemxdofexdofe), dtype=INT32)
    cdef int[:, ::1] jth = zeros((nelemxdofexdofe), dtype=INT32)
    cdef double[:, :] val = zeros((nelemxdofexdofe), dtype=FLT64)

    q0 = 0
    for i in range(nnode):
        node = inci[:, 4:].astype(INT32)
        elmlist = inci[(asarray(where(node == i + 1)))[0][:], 0]
        q1 = elmlist.size
        ith[q0 : q1 + q0] = i
        jth[q0 : q1 + q0] = elmlist - 1
        val[q0 : q1 + q0] = elmlist
        q0 = q1 + q0
    S = csc_matrix((val, (ith, jth)), shape=(nnode, nelem))
    return S

@cython.exceptval(check=False)
@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
def results_average(double[:] results_elm, int nnode, int nelem, int dofe, double[:, :] inci):
    cdef double[:, :] S = elem2nodes_conec(nnode, nelem, dofe, inci)
    cdef double[:] results_avr = zeros((nnode), dtype=FLT64)
    cdef int mm

    for mm in range(nnode):
        results_avr[mm] = mean(results_elm[(S[mm, :].nonzero())[1]])
    return results_avr

@cython.exceptval(check=False)
@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
def poly_area(double[:] a, double[:] b, double[:] c):
    cdef double result
    cdef int i
    cdef double[:] total = zeros(3, dtype=FLT64)
    cdef int N = a.shape[0]

    for i in range(N):
        prod = cross([a[i], b[i], c[i]], [a[(i + 1) % N], b[(i + 1) % N], c[(i + 1) % N]])
        total[0] += prod[0]
        total[1] += prod[1]
        total[2] += prod[2]
    result = dot(total, unit_normal(a, b, c))
    return abs(result / 2)

@cython.exceptval(check=False)
@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
def unit_normal(double[:] a, double[:] b, double[:] c):
    cdef double x, y, z, magnitude
    x = det([[1, a[1], a[2]], [1, b[1], b[2]], [1, c[1], c[2]]])
    y = det([[a[0], 1, a[2]], [b[0], 1, b[2]], [c[0], 1, c[2]]])
    z = det([[a[0], a[1], 1], [b[0], b[1], 1], [c[0], c[1], 1]])
    magnitude = (x**2 + y**2 + z**2) ** 0.5
    return (x / magnitude, y / magnitude, z / magnitude)

# @cython.exceptval(check=False)
# @cython.boundscheck(False) # turn off bounds-checking for entire function
# @cython.wraparound(False)  # turn off negative index wrapping for entire function
# def search_edgex(double edge_coord, double[:, :] coord, double erro):
#     cdef double [:] dif = abs(edge_coord*ones_like(coord[:, 1]) - coord[:, 1])
#     cdef double [:] erroa = erro*ones_like(dif)
#     cdef bint [:] idx = less(dif, erroa) #array(where(dif < erroa)[0][:]).astype(INT32) #(asarray(where(dif < erroa)))[0][:]
#     cdef int [:] node = coord[idx, 0].astype(INT32)
#     return node

# @cython.exceptval(check=False)
# @cython.boundscheck(False) # turn off bounds-checking for entire function
# @cython.wraparound(False)  # turn off negative index wrapping for entire function
# def search_edgey(double edge_coord, double[:, :] coord, double erro):
#     cdef double dif = abs(edge_coord * ones_like(coord[:, 2]) - coord[:, 2])
#     cdef int node_posi = (asarray(where(dif < erro)))[0][:]
#     cdef int[:] node = coord[node_posi, 0].astype(INT32)
#     return node

# @cython.exceptval(check=False)
# @cython.boundscheck(False) # turn off bounds-checking for entire function
# @cython.wraparound(False)  # turn off negative index wrapping for entire function
# def search_edgez(double edge_coord, double[:, :] coord, double erro):
#     cdef double dif = abs(edge_coord * ones_like(coord[:, 3]) - coord[:, 3])
#     cdef int node_posi = (asarray(where(dif < erro)))[0][:]
#     cdef int[:] node = coord[node_posi, 0].astype(INT32)
#     return node

# @cython.exceptval(check=False)
# @cython.boundscheck(False) # turn off bounds-checking for entire function
# @cython.wraparound(False)  # turn off negative index wrapping for entire function
# def search_surfxy(double[:] orthg_coord, double[:, :] coord, double erro):
#     cdef double dif = abs(orthg_coord * ones_like(coord[:, 3]) - coord[:, 3])
#     cdef int node_posi = (asarray(where(dif < erro)))[0][:]
#     cdef int[:] node = coord[node_posi, 0].astype(INT32)
#     return node

# @cython.exceptval(check=False)
# @cython.boundscheck(False) # turn off bounds-checking for entire function
# @cython.wraparound(False)  # turn off negative index wrapping for entire function
# def search_surfyz(double[:] orthg_coord, double[:, :] coord, double erro):
#     cdef double dif = abs(orthg_coord * ones_like(coord[:, 1]) - coord[:, 1])
#     cdef int node_posi = (asarray(where(dif < erro)))[0][:]
#     cdef int[:] node = coord[node_posi, 0].astype(INT32)
#     return node

# @cython.exceptval(check=False)
# @cython.boundscheck(False) # turn off bounds-checking for entire function
# @cython.wraparound(False)  # turn off negative index wrapping for entire function
# def search_surfzx(double[:] orthg_coord, double[:, :] coord, double erro):
#     cdef double dif = abs(orthg_coord * ones_like(coord[:, 2]) - coord[:, 2])
#     cdef int node_posi = (asarray(where(dif < erro)))[0][:]
#     cdef int[:] node = coord[node_posi, 0].astype(INT32)
#     return node

# @cython.exceptval(check=False)
# @cython.boundscheck(False) # turn off bounds-checking for entire function
# @cython.wraparound(False)  # turn off negative index wrapping for entire function
# def search_nodexyz(double[:] node_coordX, double[:] node_coordY, double[:] node_coordZ, double[:, :] coord, double erro):
#     cdef double difx = abs(node_coordX * ones_like(coord[:, 1]) - coord[:, 1])
#     cdef double dify = abs(node_coordY * ones_like(coord[:, 2]) - coord[:, 2])
#     cdef double difz = abs(node_coordZ * ones_like(coord[:, 3]) - coord[:, 3])
#     cdef int node_posi = (asarray(where((difx < erro) & (dify < erro) & (difz < erro))))[0][:]
#     cdef int[:] node = coord[node_posi, 0].astype(INT32)
#     return node

# @cython.exceptval(check=False)
# @cython.boundscheck(False) # turn off bounds-checking for entire function
# @cython.wraparound(False)  # turn off negative index wrapping for entire function
# def nodes_from_regions(dict regionlist):
#     cdef list points=[], edges=[], surfs=[], regions=[]
#     cdef dict pointlist, edgelist, surflist
#     cdef int pp, ee, ss
#     cdef int[:] nodes

#     pointlist = regionlist["point"]
#     for pp in range(len(pointlist)):
#         nodes = unique(array(pointlist[str(pp + 1)]["nodes"]))
#         points.append([pp + 1, nodes])

#     edgelist = regionlist["line"]
#     for ee in range(len(edgelist)):
#         nodes = unique(array(edgelist[str(ee + 1)]["nodes"]))
#         edges.append([ee + 1, nodes])

#     surflist = regionlist["plane"]
#     for ss in range(len(surflist)):
#         nodes = unique(array(surflist[str(ss + 1)]["nodes"]))
#         surfs.append([ss + 1, nodes])

#     regions.append(["POINT_", points])
#     regions.append(["LINE_", edges])
#     regions.append(["PLANE_", surfs])

#     return regions


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

cpdef __gauss_points_quad4(int npp):
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
              [0.7745966692414834, 0.7745966692414834]]
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
        