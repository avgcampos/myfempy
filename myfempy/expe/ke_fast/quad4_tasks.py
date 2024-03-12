import numpy as np
INT32 = np.int32
FLT64 = np.float64

                    
def ShapeFunctions(r_coord, nodedof):
    shape_function =  __N(r_coord)
    mat_N = np.zeros((nodedof, 4*nodedof), dtype=FLT64) 
    for block in range(4):
        for dof in range(nodedof):
            mat_N[dof, block*nodedof+dof] = shape_function[0, block]
    return  mat_N


def DiffShapeFuntion(r_coord, nodedof):
    diff_shape_function = __diffN(r_coord)
    mat_diff_N = np.zeros((2*nodedof, 4*nodedof), dtype=FLT64) 
    for block in range(4):
        for dof in range(nodedof):
            mat_diff_N[nodedof*dof-dof*(nodedof-2), block*nodedof+dof] = diff_shape_function[0, block]
            mat_diff_N[nodedof*dof-dof*(nodedof-2)+1, block*nodedof+dof] = diff_shape_function[1, block]
    return mat_diff_N
    
    
def Jacobian(r_coord, element_coord):  
    diffN = __diffN(r_coord)  
    jac = np.dot(diffN, element_coord)
    return jac
    

def invJacobi(r_coord, element_coord, nodedof):
    J = Jacobian(r_coord, element_coord)
    invJ = __inverse(J.flatten())
    mat_invJ = np.zeros((2*nodedof, 2*nodedof), dtype=FLT64)  
    mat_invJ[0, 0] = invJ[0, 0]
    mat_invJ[0, 1] = invJ[0, 1]
    mat_invJ[1, 0] = invJ[1, 0]
    mat_invJ[1, 1] = invJ[1, 1]
    mat_invJ[2, 2] = invJ[0, 0]
    mat_invJ[2, 3] = invJ[0, 1]
    mat_invJ[3, 2] = invJ[1, 0]
    mat_invJ[3, 3] = invJ[1, 1]     
    return mat_invJ


def detJacobi(r_coord, element_coord):
    J = Jacobian(r_coord, element_coord)
    detJ = __determinant(J.flatten())
    return detJ


def NodeList(inci, element_number):
    noi = int(inci[element_number, 4])
    noj = int(inci[element_number, 5])
    nok = int(inci[element_number, 6])
    nol = int(inci[element_number, 7])
    node_list = [noi, noj, nok, nol]                  
    return node_list
            
def NodeCoord(coord, node_list):
    noi = node_list[0]
    noj = node_list[1]
    nok = node_list[2]
    nol = node_list[3]
    xi = coord[noi - 1, 1]
    yi = coord[noi - 1, 2]
    xj = coord[noj - 1, 1]
    yj = coord[noj - 1, 2]
    xk = coord[nok - 1, 1]
    yk = coord[nok - 1, 2]
    xl = coord[nol - 1, 1]
    yl = coord[nol - 1, 2]
    element_coord = np.array([[xi, yi], [xj, yj], [xk, yk], [xl, yl]], dtype=FLT64)
    return element_coord


def LocKey(node_list, nodedof):
    shape_key = np.zeros(4*nodedof, dtype=INT32)
    for node in range(len(node_list)):
        for dof in range(nodedof):
            shape_key[nodedof*node+dof] = nodedof * node_list[node] - (nodedof-dof)
    return shape_key


def __N(r_coord):
    N = np.zeros((1,4), dtype=FLT64)
    N[0, 0] = 0.25*(1-r_coord[0])*(1-r_coord[1])
    N[0, 1] = 0.25*(1+r_coord[0])*(1-r_coord[1])
    N[0, 2] = 0.25*(1+r_coord[0])*(1+r_coord[1])
    N[0, 3] = 0.25*(1-r_coord[0])*(1+r_coord[1])
    return N

# @profile
def __diffN(r):
    r1 = r[1] 
    r0 = r[0]
    dN = np.zeros((2,4), dtype=FLT64)
    dN[0, 0] = 0.25*(-1.0+r1)
    dN[0, 1] = 0.25*(1.0-r1)
    dN[0, 2] = 0.25*(1.0+r1)
    dN[0, 3] = 0.25*(-1.0-r1)
    dN[1, 0] = 0.25*(-1.0+r0)
    dN[1, 1] = 0.25*(-1.0-r0)
    dN[1, 2] = 0.25*(1.0+r0)
    dN[1, 3] = 0.25*(1.0-r0)
    return dN

def __determinant(A):
    detA = A[0]*A[3]-A[1]*A[2]
    return detA


def __inverse(A):
    invA = 1/(A[0]*A[3]-A[1]*A[2])*np.array([[A[3], -A[1]], [-A[2], A[0]]])
    return invA