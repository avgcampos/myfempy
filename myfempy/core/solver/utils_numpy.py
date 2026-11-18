import numpy as np
import scipy as sp


def getVectorization(ith, jth, val, loc, matrix, ee, elemdof):

    ith[(elemdof * elemdof) * ee : (elemdof * elemdof) * (ee + 1)] = np.tile(
        loc.reshape(1, elemdof).T, (1, elemdof)
    ).flatten("F")
    jth[(elemdof * elemdof) * ee : (elemdof * elemdof) * (ee + 1)] = np.transpose(
        np.tile(loc.reshape(1, elemdof).T, (1, elemdof))
    ).flatten("F")
    val[(elemdof * elemdof) * ee : (elemdof * elemdof) * (ee + 1)] = matrix.flatten("F")
    return ith, jth, val

def getLoadAssembler(loadaply, nodetot, nodedof):
    """
    getLoadAssembler Assembler module <ConcreteClassService>

    Returns:
        vector sparse
    """
    unique_steps = np.unique(loadaply[:, 3])
    steps = max(1, len(unique_steps))

    forcevec = np.zeros((nodedof * nodetot, steps))

    for fstep in range(len(unique_steps)):
        forceaply = loadaply[loadaply[:, 3] == fstep + 1]
        
        for ii in range(forceaply.shape[0]):
            node = int(forceaply[ii, 0])
            dof = int(forceaply[ii, 1])
            val = forceaply[ii, 2]

            # Garante que o DOF informado está dentro do limite do elemento
            if 1 <= dof <= nodedof:
                gdlload = (node - 1) * nodedof + (dof - 1)
                forcevec[gdlload, fstep] += val

    return forcevec

import numpy as np

# Dirichlet Homogeneous https://en.wikipedia.org/wiki/Dirichlet_boundary_condition
def getConstrains(constrains, nodetot, nodedof):
    """
    getConstrains Constrain module <ConcreteClassService>

    Returns:
        freedof, fixedof, constdof (arrays 1D com índices em base 0)
    """
    total_dof = nodedof * nodetot
    fixed_mask = np.zeros(total_dof, dtype=bool)
    const_mask = np.zeros(total_dof, dtype=bool)

    for ii in range(len(constrains)):
        no = int(constrains[ii, 0])
        dof = int(constrains[ii, 1])
        val = constrains[ii, 2]

        # Se dof == 0, restringe todos os DOFs do nó; senão, apenas o DOF específico
        if dof == 0:
            dofs_to_constrain = range(1, nodedof + 1)
        elif 1 <= dof <= nodedof:
            dofs_to_constrain = (dof,)
        else:
            continue  # Ignora DOFs inválidos para o nodedof atual

        for d in dofs_to_constrain:
            idx = (no - 1) * nodedof + (d - 1)
            if val == 0.0:
                fixed_mask[idx] = True
            else:
                const_mask[idx] = True

    alldof = np.arange(total_dof, dtype=int)
    fixedof = alldof[fixed_mask]
    constdof = alldof[const_mask]
    freedof = alldof[~(fixed_mask | const_mask)]

    return freedof, fixedof, constdof

import numpy as np

# Dirichlet Non-Homogeneous
def getDirichletNH(constrains, nodetot, nodedof):
    """
    getDirichletNH Assembler module <ConcreteClassService>

    Returns:
        Uc (ndarray): vetor/matriz de deslocamentos prescritos
    """
    unique_steps = np.unique(constrains[:, 3][constrains[:, 3] != 0])
    steps = max(1, len(unique_steps))

    Uc = np.zeros((nodedof * nodetot, steps), dtype=np.float64)  # solution constrains

    for cstep in range(len(unique_steps)):
        forceaply = constrains[constrains[:, 3] == cstep + 1]

        for ii in range(forceaply.shape[0]):
            node = int(forceaply[ii, 0])
            dof = int(forceaply[ii, 1])
            val = forceaply[ii, 2]

            if 1 <= dof <= nodedof:
                gdlload = (node - 1) * nodedof + (dof - 1)
                Uc[gdlload, cstep] = val

    return Uc

# https://en.wikipedia.org/wiki/Rotation_matrix
def getRotationMatrix(node_list, coord, ndof, dim=None, origin=None, ref_vec=None):
    """
    Constrói a matriz de rotação global em blocos esparsos para sistemas 2D ou 3D.

    Parâmetros:
        node_list : array 1D com os IDs dos nós (base 1).
        coord     : array com coordenadas [ID, X, Y] (2D) ou [ID, X, Y, Z] (3D).
        ndof      : número total de graus de liberdade (nodedof * nodetot).
        dim       : 2 ou 3. Se None, infere a partir das colunas de `coord`.
        origin    : coordenadas da origem [x0, y0] ou [x0, y0, z0]. Padrão: zeros.
        ref_vec   : vetor auxiliar 3D para definir o plano de referência. Padrão: [0, 0, 1].
    """
    # Infere a dimensão automaticamente se não for informada (coord tem coluna ID + 2 ou 3 coords)
    if dim is None:
        dim = 3 if coord.shape[1] >= 4 else 2

    if origin is None:
        origin = np.zeros(dim, dtype=float)
    else:
        origin = np.asarray(origin, dtype=float)

    if ref_vec is None:
        ref_vec = np.array([0.0, 0.0, 1.0])
    else:
        ref_vec = np.asarray(ref_vec, dtype=float)

    n_nodes = node_list.shape[0]
    nodedof = ndof // n_nodes  # Suporta 2 (2D), 3 (3D) ou 6 (3D com rotações)

    RM = sp.sparse.lil_matrix((ndof, ndof), dtype=np.float64)

    for n in range(n_nodes):
        nol = int(node_list[n]) - 1
        pos = coord[nol, 1 : 1 + dim] - origin
        norm = np.linalg.norm(pos)

        # Se o nó está na origem (norm == 0), mantém a matriz identidade (sem rotação)
        if np.isclose(norm, 0.0):
            R_node = np.eye(dim)
        elif dim == 2:
            c_the = pos[0] / norm
            s_the = pos[1] / norm
            R_node = np.array([
                [c_the, -s_the],
                [s_the,  c_the]
            ])
        elif dim == 3:
            e1 = pos / norm
            # Evita produto vetorial nulo se e1 for paralelo a ref_vec
            v_aux = ref_vec if abs(np.dot(e1, ref_vec)) < 0.999 else np.array([0.0, 1.0, 0.0])
            e2 = np.cross(v_aux, e1)
            e2 /= np.linalg.norm(e2)
            e3 = np.cross(e1, e2)
            R_node = np.column_stack((e1, e2, e3))
        else:
            raise ValueError("A dimensão `dim` deve ser 2 ou 3.")

        # Preenche o(s) bloco(s) na diagonal da matriz esparsa
        base = nodedof * n
        if nodedof == 6 and dim == 3:
            # Para elementos de pórtico/casca 3D (3 translações + 3 rotações)
            RM[base : base + 3, base : base + 3] = R_node
            RM[base + 3 : base + 6, base + 3 : base + 6] = R_node
        else:
            RM[base : base + dim, base : base + dim] = R_node

    return RM.tocsr()

# def getLoadAssembler(loadaply, nodetot, nodedof):
#     """
#     getLoadAssembler Assembler module <ConcreteClassService>

#     Returns:
#         vector sparse
#     """

#     steps = len(np.unique(loadaply[:, 3]))
#     if steps == 0:
#         steps = 1
#     else:
#         pass

#     forcevec = np.zeros((nodedof * nodetot, steps))

#     for fstep in range(len(np.unique(loadaply[:, 3]))):
#         forceaply = loadaply[np.where(loadaply[:, 3] == fstep + 1), :][0]
#         nload = forceaply.shape[0]

#         if nodedof == 1:
#             for ii in range(nload):
#                 if int(forceaply[ii, 1]) == 1:
#                     gdlload = int(nodedof * forceaply[ii, 0] - (nodedof))
#                     forcevec[gdlload, fstep] += forceaply[ii, 2]
#                 elif int(forceaply[ii, 1]) == 2:
#                     gdlload = int(nodedof * forceaply[ii, 0] - (nodedof))
#                     forcevec[gdlload, fstep] += forceaply[ii, 2]

#         elif nodedof == 2:
#             for ii in range(nload):
#                 if int(forceaply[ii, 1]) == 1:
#                     gdlload = int(nodedof * forceaply[ii, 0] - (nodedof))
#                     forcevec[gdlload, fstep] += forceaply[ii, 2]

#                 elif int(forceaply[ii, 1]) == 2:
#                     gdlload = int(nodedof * forceaply[ii, 0] - (nodedof - 1))
#                     forcevec[gdlload, fstep] += forceaply[ii, 2]

#         elif nodedof == 3:
#             for ii in range(nload):
#                 if int(forceaply[ii, 1]) == 1:
#                     gdlload = int(nodedof * forceaply[ii, 0] - (nodedof))
#                     forcevec[gdlload, fstep] += forceaply[ii, 2]

#                 elif int(forceaply[ii, 1]) == 2:
#                     gdlload = int(nodedof * forceaply[ii, 0] - (nodedof - 1))
#                     forcevec[gdlload, fstep] += forceaply[ii, 2]

#                 elif int(forceaply[ii, 1]) == 3:
#                     gdlload = int(nodedof * forceaply[ii, 0] - (nodedof - 2))
#                     forcevec[gdlload, fstep] += forceaply[ii, 2]

#         elif nodedof == 6:
#             for ii in range(nload):
#                 if int(forceaply[ii, 1]) == 1:
#                     gdlload = int(nodedof * forceaply[ii, 0] - (nodedof))
#                     forcevec[gdlload, fstep] += forceaply[ii, 2]

#                 elif int(forceaply[ii, 1]) == 2:
#                     gdlload = int(nodedof * forceaply[ii, 0] - (nodedof - 1))
#                     forcevec[gdlload, fstep] += forceaply[ii, 2]

#                 elif int(forceaply[ii, 1]) == 3:
#                     gdlload = int(nodedof * forceaply[ii, 0] - (nodedof - 2))
#                     forcevec[gdlload, fstep] += forceaply[ii, 2]

#                 elif int(forceaply[ii, 1]) == 4:
#                     gdlload = int(nodedof * forceaply[ii, 0] - (nodedof - 3))
#                     forcevec[gdlload, fstep] += forceaply[ii, 2]

#                 elif int(forceaply[ii, 1]) == 5:
#                     gdlload = int(nodedof * forceaply[ii, 0] - (nodedof - 4))
#                     forcevec[gdlload, fstep] += forceaply[ii, 2]

#                 elif int(forceaply[ii, 1]) == 6:
#                     gdlload = int(nodedof * forceaply[ii, 0] - (nodedof - 5))
#                     forcevec[gdlload, fstep] += forceaply[ii, 2]
#         else:
#             pass
#     # forcevec = sparse.csc_matrix(forcevec)
#     return forcevec


# # Dirichlet Homogeneous https://en.wikipedia.org/wiki/Dirichlet_boundary_condition
# def getConstrains(constrains, nodetot, nodedof):
#     """
#     getConstrains Constrain module <ConcreteClassService>

#     Returns:
#         _description_
#     """
#     ntbc = len(constrains)
#     fixedof = np.zeros((nodedof * nodetot, 1), dtype=int)
#     constdof = np.zeros((nodedof * nodetot, 1), dtype=int)

#     if nodedof == 1:
#         for ii in range(ntbc):
#             no = int(constrains[ii, 0])
#             if int(constrains[ii, 1]) == 1:
#                 if constrains[ii, 2] == 0.0:
#                     fixedof[nodedof * no - 1, 0] = nodedof * no
#                 else:
#                     constdof[nodedof * no - 1, 0] = nodedof * no
#             else:
#                 if constrains[ii, 2] == 0.0:
#                     fixedof[nodedof * no - 1, 0] = nodedof * no
#                 else:
#                     constdof[nodedof * no - 1, 0] = nodedof * no

#     elif nodedof == 2:
#         for ii in range(ntbc):
#             no = int(constrains[ii, 0])
#             if int(constrains[ii, 1]) == 1:
#                 if constrains[ii, 2] == 0.0:
#                     fixedof[nodedof * no - 2, 0] = nodedof * no - 1
#                 else:
#                     constdof[nodedof * no - 2, 0] = nodedof * no - 1
#             elif int(constrains[ii, 1]) == 2:
#                 if constrains[ii, 2] == 0.0:
#                     fixedof[nodedof * no - 1, 0] = nodedof * no
#                 else:
#                     constdof[nodedof * no - 1, 0] = nodedof * no
#             else:  # int(constrains[ii, 1]) == 0:
#                 if constrains[ii, 2] == 0.0:
#                     fixedof[nodedof * no - 2, 0] = nodedof * no - 1
#                     fixedof[nodedof * no - 1, 0] = nodedof * no
#                 else:
#                     constdof[0, nodedof * no - 2, 0] = nodedof * no - 1
#                     constdof[0, nodedof * no - 1, 0] = nodedof * no

#     elif nodedof == 3:
#         for ii in range(ntbc):
#             no = int(constrains[ii, 0])
#             if int(constrains[ii, 1]) == 1:
#                 if constrains[ii, 2] == 0.0:
#                     fixedof[nodedof * no - 3, 0] = nodedof * no - 2
#                 else:
#                     constdof[nodedof * no - 3, 0] = nodedof * no - 2
#             elif int(constrains[ii, 1]) == 2:
#                 if constrains[ii, 2] == 0.0:
#                     fixedof[nodedof * no - 2, 0] = nodedof * no - 1
#                 else:
#                     constdof[nodedof * no - 2, 0] = nodedof * no - 1
#             elif int(constrains[ii, 1]) == 3:
#                 if constrains[ii, 2] == 0.0:
#                     fixedof[nodedof * no - 1, 0] = nodedof * no
#                 else:
#                     constdof[nodedof * no - 1, 0] = nodedof * no
#             else:  # int(constrains[ii, 1]) == 0:
#                 if constrains[ii, 2] == 0.0:
#                     fixedof[nodedof * no - 3, 0] = nodedof * no - 2
#                     fixedof[nodedof * no - 2, 0] = nodedof * no - 1
#                     fixedof[nodedof * no - 1, 0] = nodedof * no
#                 else:
#                     constdof[nodedof * no - 3, 0] = nodedof * no - 2
#                     constdof[nodedof * no - 2, 0] = nodedof * no - 1
#                     constdof[nodedof * no - 1, 0] = nodedof * no

#     elif nodedof == 6:
#         for ii in range(ntbc):
#             no = int(constrains[ii, 0])
#             if int(constrains[ii, 1]) == 1:
#                 if constrains[ii, 2] == 0.0:
#                     fixedof[nodedof * no - 6, 0] = nodedof * no - 5
#                 else:
#                     constdof[nodedof * no - 6, 0] = nodedof * no - 5

#             elif int(constrains[ii, 1]) == 2:
#                 if constrains[ii, 2] == 0.0:
#                     fixedof[nodedof * no - 5, 0] = nodedof * no - 4
#                 else:
#                     constdof[nodedof * no - 5, 0] = nodedof * no - 4

#             elif int(constrains[ii, 1]) == 3:
#                 if constrains[ii, 2] == 0.0:
#                     fixedof[nodedof * no - 4, 0] = nodedof * no - 3
#                 else:
#                     constdof[nodedof * no - 4, 0] = nodedof * no - 3

#             elif int(constrains[ii, 1]) == 4:
#                 if constrains[ii, 2] == 0.0:
#                     fixedof[nodedof * no - 3, 0] = nodedof * no - 2
#                 else:
#                     constdof[nodedof * no - 3, 0] = nodedof * no - 2

#             elif int(constrains[ii, 1]) == 5:
#                 if constrains[ii, 2] == 0.0:
#                     fixedof[nodedof * no - 2, 0] = nodedof * no - 1
#                 else:
#                     constdof[nodedof * no - 2, 0] = nodedof * no - 1

#             elif int(constrains[ii, 1]) == 6:
#                 if constrains[ii, 2] == 0.0:
#                     fixedof[nodedof * no - 1, 0] = nodedof * no
#                 else:
#                     constdof[nodedof * no - 1, 0] = nodedof * no

#             else:  # int(constrains[ii, 1]) == 0:
#                 if constrains[ii, 2] == 0.0:
#                     fixedof[nodedof * no - 6, 0] = nodedof * no - 5
#                     fixedof[nodedof * no - 5, 0] = nodedof * no - 4
#                     fixedof[nodedof * no - 4, 0] = nodedof * no - 3
#                     fixedof[nodedof * no - 3, 0] = nodedof * no - 2
#                     fixedof[nodedof * no - 2, 0] = nodedof * no - 1
#                     fixedof[nodedof * no - 1, 0] = nodedof * no
#                 else:
#                     constdof[nodedof * no - 6, 0] = nodedof * no - 5
#                     constdof[nodedof * no - 5, 0] = nodedof * no - 4
#                     constdof[nodedof * no - 4, 0] = nodedof * no - 3
#                     constdof[nodedof * no - 3, 0] = nodedof * no - 2
#                     constdof[nodedof * no - 2, 0] = nodedof * no - 1
#                     constdof[nodedof * no - 1, 0] = nodedof * no

#     fixedof = fixedof[np.nonzero(fixedof)]
#     fixedof = fixedof - np.ones_like(fixedof)
#     constdof = constdof[np.nonzero(constdof)]
#     constdof = constdof - np.ones_like(constdof)
#     alldof = np.arange(0, nodedof * nodetot, 1, int)
#     freedof = np.setdiff1d(alldof, fixedof)
#     freedof = np.setdiff1d(freedof, constdof)
#     return freedof, fixedof, constdof


# # Dirichlet Non-Homogeneous
# def getDirichletNH(constrains, nodetot, nodedof):
#     """
#     getLoadAssembler Assembler module <ConcreteClassService>

#     Returns:
#         vector sparse
#     """

#     steps = len(np.unique(constrains[:, 3][constrains[:, 3] != 0]))
#     if steps == 0:
#         steps = 1
#     else:
#         pass

#     Uc = np.zeros((nodedof * nodetot, steps), dtype=np.float64)  # solution constrains

#     for cstep in range(len(np.unique(constrains[:, 3][constrains[:, 3] != 0]))):
#         forceaply = constrains[np.where(constrains[:, 3] == cstep + 1), :][0]
#         nload = forceaply.shape[0]

#         if nodedof == 1:
#             for ii in range(nload):
#                 if int(forceaply[ii, 1]) == 1:
#                     gdlload = int(nodedof * forceaply[ii, 0] - (nodedof))
#                     Uc[gdlload, cstep] = forceaply[ii, 2]

#         elif nodedof == 2:
#             for ii in range(nload):
#                 if int(forceaply[ii, 1]) == 1:
#                     gdlload = int(nodedof * forceaply[ii, 0] - (nodedof))
#                     Uc[gdlload, cstep] = forceaply[ii, 2]

#                 elif int(forceaply[ii, 1]) == 2:
#                     gdlload = int(nodedof * forceaply[ii, 0] - (nodedof - 1))
#                     Uc[gdlload, cstep] = forceaply[ii, 2]

#         elif nodedof == 3:
#             for ii in range(nload):
#                 if int(forceaply[ii, 1]) == 1:
#                     gdlload = int(nodedof * forceaply[ii, 0] - (nodedof))
#                     Uc[gdlload, cstep] = forceaply[ii, 2]

#                 elif int(forceaply[ii, 1]) == 2:
#                     gdlload = int(nodedof * forceaply[ii, 0] - (nodedof - 1))
#                     Uc[gdlload, cstep] = forceaply[ii, 2]

#                 elif int(forceaply[ii, 1]) == 3:
#                     gdlload = int(nodedof * forceaply[ii, 0] - (nodedof - 2))
#                     Uc[gdlload, cstep] = forceaply[ii, 2]

#         elif nodedof == 6:
#             for ii in range(nload):
#                 if int(forceaply[ii, 1]) == 1:
#                     gdlload = int(nodedof * forceaply[ii, 0] - (nodedof))
#                     Uc[gdlload, cstep] = forceaply[ii, 2]

#                 elif int(forceaply[ii, 1]) == 2:
#                     gdlload = int(nodedof * forceaply[ii, 0] - (nodedof - 1))
#                     Uc[gdlload, cstep] = forceaply[ii, 2]

#                 elif int(forceaply[ii, 1]) == 3:
#                     gdlload = int(nodedof * forceaply[ii, 0] - (nodedof - 2))
#                     Uc[gdlload, cstep] = forceaply[ii, 2]

#                 elif int(forceaply[ii, 1]) == 4:
#                     gdlload = int(nodedof * forceaply[ii, 0] - (nodedof - 3))
#                     Uc[gdlload, cstep] = forceaply[ii, 2]

#                 elif int(forceaply[ii, 1]) == 5:
#                     gdlload = int(nodedof * forceaply[ii, 0] - (nodedof - 4))
#                     Uc[gdlload, cstep] = forceaply[ii, 2]

#                 elif int(forceaply[ii, 1]) == 6:
#                     gdlload = int(nodedof * forceaply[ii, 0] - (nodedof - 5))
#                     Uc[gdlload, cstep] = forceaply[ii, 2]

#         else:
#             pass
#     return Uc


# https://en.wikipedia.org/wiki/Rotation_matrix
def getRotationMatrix2D(node_list, coord, ndof):
    # Initialize RM as a sparse matrix
    RM = sp.sparse.lil_matrix((ndof, ndof))
    for n in range(node_list.shape[0]):
        nol = int(node_list[n] - 1)
        RonX = coord[nol, 1]  # - Og[0]
        RonY = coord[nol, 2]  # - Og[1]
        Ron = np.sqrt(RonX**2 + RonY**2)
        S_the = RonY / Ron
        C_the = RonX / Ron
        # Assign values to the sparse matrix
        RM[2 * n, 2 * n] = C_the
        RM[2 * n, 2 * n + 1] = -S_the
        RM[2 * n + 1, 2 * n] = S_the
        RM[2 * n + 1, 2 * n + 1] = C_the

    # Convert to CSR format for more efficient arithmetic and matrix-vector operations
    return RM.tocsr()

# def getRotationMatrix3D(node_list, coord, ndof):
#     """
#     https://en.wikipedia.org/wiki/Rotation_matrix

#     Arguments:
#         node_list -- _description_
#         coord -- _description_
#         ndof -- _description_

#     Returns:
#         RM
#     """
#     # Inicialize a RM como uma matriz esparsa
#     RM = lil_matrix((ndof, ndof))
#     for n in range(node_list.shape[0]):
#         nol = int(node_list[n] - 1)
#         RonX = coord[nol, 1]
#         RonY = coord[nol, 2]
#         RonZ = coord[nol, 3]
#         Ron = np.sqrt(RonX**2 + RonY**2 + RonZ**2)
#         S_phi = RonY / Ron
#         C_phi = RonX / Ron
#         S_theta = RonZ / Ron
#         C_theta = np.sqrt(RonX**2 + RonY**2) / Ron

#         # Atribuir valores à matriz esparsa
#         RM[3 * n, 3 * n] = C_theta
#         RM[3 * n, 3 * n + 1] = -S_phi * S_theta
#         RM[3 * n, 3 * n + 2] = C_phi * S_theta
#         RM[3 * n + 1, 3 * n] = S_phi
#         RM[3 * n + 1, 3 * n + 1] = C_phi * C_theta
#         RM[3 * n + 1, 3 * n + 2] = -S_theta
#         RM[3 * n + 2, 3 * n] = -S_theta
#         RM[3 * n + 2, 3 * n + 1] = C_theta
#         RM[3 * n + 2, 3 * n + 2] = C_theta * C_phi

#     # Converter para o formato CSR para operações aritméticas e de matriz-vetor mais eficientes
#     return RM.tocsr()