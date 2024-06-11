from __future__ import annotations

import numpy as np
from scipy import sparse

# from myfempy.core.alglin import linsolve_spsolve
from myfempy.core.solver.assembler import Assembler, getLoc, getMatrix
from myfempy.core.solver.assemblerfull_numpy_v1 import getMatrixAssembler_Full

# from myfempy.core.solver.assemblerfull_cython_v4 import getMatrixAssembler_Full


class AssemblerFULL(Assembler):

    """
    Assembler Full System Class <ConcreteClassService>
    """

    def getMatrixAssembler(
        Model, inci, coord, tabmat, tabgeo, intgauss, type_assembler, MP
    ):
        """
        getMatrixAssembler Assembler module <ConcreteClassService>

        Returns:
            matrix sparse
        """
        elem_set = Model.element.getElementSet()
        nodedof = len(elem_set["dofs"]["d"])
        shape_set = Model.shape.getShapeSet()
        nodecon = len(shape_set["nodes"])
        elemdof = nodecon * nodedof
        nodetot = coord.shape[0]
        sdof = nodedof * nodetot

        # if MP>0:
        #     ith, jth, val = getMatrixAssemblerSym_cy_mp(Model, inci, coord, tabmat, tabgeo, elemdof,  intgauss, type_assembler, MP)
        # else:
        ith, jth, val = getMatrixAssembler_Full(
            Model, inci, coord, tabmat, tabgeo, elemdof, intgauss, type_assembler
        )

        A_sp_scipy_csc = sparse.csc_matrix(
            (val, (ith, jth)),
            shape=(sdof, sdof),
        )
        return A_sp_scipy_csc

    def getLoadAssembler(loadaply, nodetot, nodedof):
        """
        getLoadAssembler Assembler module <ConcreteClassService>

        Returns:
            vector sparse
        """

        steps = len(np.unique(loadaply[:, 3]))
        if steps == 0:
            steps = 1
        else:
            pass

        forcevec = np.zeros((nodedof * nodetot, steps))

        for fstep in range(len(np.unique(loadaply[:, 3]))):
            forceaply = loadaply[np.where(loadaply[:, 3] == fstep + 1), :][0]
            nload = forceaply.shape[0]

            if nodedof == 1:
                for ii in range(nload):
                    if int(forceaply[ii, 1]) == 1:
                        gdlload = int(nodedof * forceaply[ii, 0] - (nodedof - 1))
                        forcevec[gdlload, fstep] += forceaply[ii, 2]

            elif nodedof == 2:
                for ii in range(nload):
                    if int(forceaply[ii, 1]) == 1:
                        gdlload = int(nodedof * forceaply[ii, 0] - (nodedof))
                        forcevec[gdlload, fstep] += forceaply[ii, 2]

                    elif int(forceaply[ii, 1]) == 2:
                        gdlload = int(nodedof * forceaply[ii, 0] - (nodedof - 1))
                        forcevec[gdlload, fstep] += forceaply[ii, 2]

            elif nodedof == 3:
                for ii in range(nload):
                    if int(forceaply[ii, 1]) == 1:
                        gdlload = int(nodedof * forceaply[ii, 0] - (nodedof))
                        forcevec[gdlload, fstep] += forceaply[ii, 2]

                    elif int(forceaply[ii, 1]) == 2:
                        gdlload = int(nodedof * forceaply[ii, 0] - (nodedof - 1))
                        forcevec[gdlload, fstep] += forceaply[ii, 2]

                    elif int(forceaply[ii, 1]) == 3:
                        gdlload = int(nodedof * forceaply[ii, 0] - (nodedof - 2))
                        forcevec[gdlload, fstep] += forceaply[ii, 2]

            else:
                pass
        # forcevec = sparse.csc_matrix(forcevec)
        return forcevec

    # Dirichlet Homogeneous https://en.wikipedia.org/wiki/Dirichlet_boundary_condition
    def getConstrains(constrains, nodetot, nodedof):
        """
        getConstrains Constrain module <ConcreteClassService>

        Returns:
            _description_
        """

        ntbc = len(constrains)
        fixedof = np.zeros((nodedof * nodetot, 1), dtype=int)
        constdof = np.zeros((nodedof * nodetot, 1), dtype=int)

        if nodedof == 1:
            for ii in range(ntbc):
                no = int(constrains[ii, 0])
                if int(constrains[ii, 1]) == 1:
                    if constrains[ii, 2] == 0.0:
                        fixedof[nodedof * no - 1, 0] = nodedof * no
                    else:
                        constdof[nodedof * no - 1, 0] = nodedof * no

        elif nodedof == 2:
            for ii in range(ntbc):
                no = int(constrains[ii, 0])
                if int(constrains[ii, 1]) == 0:
                    if constrains[ii, 2] == 0.0:
                        fixedof[nodedof * no - 2, 0] = nodedof * no - 1
                        fixedof[nodedof * no - 1, 0] = nodedof * no
                    else:
                        constdof[0, nodedof * no - 2, 0] = nodedof * no - 1
                        constdof[0, nodedof * no - 1, 0] = nodedof * no
                elif int(constrains[ii, 1]) == 1:
                    if constrains[ii, 2] == 0.0:
                        fixedof[nodedof * no - 2, 0] = nodedof * no - 1
                    else:
                        constdof[nodedof * no - 2, 0] = nodedof * no - 1
                elif int(constrains[ii, 1]) == 2:
                    if constrains[ii, 2] == 0.0:
                        fixedof[nodedof * no - 1, 0] = nodedof * no
                    else:
                        constdof[nodedof * no - 1, 0] = nodedof * no

        elif nodedof == 3:
            for ii in range(ntbc):
                no = int(constrains[ii, 0])
                if int(constrains[ii, 1]) == 0:
                    if constrains[ii, 2] == 0.0:
                        fixedof[nodedof * no - 3, 0] = nodedof * no - 2
                        fixedof[nodedof * no - 2, 0] = nodedof * no - 1
                        fixedof[nodedof * no - 1, 0] = nodedof * no
                    else:
                        constdof[nodedof * no - 3, 0] = nodedof * no - 2
                        constdof[nodedof * no - 2, 0] = nodedof * no - 1
                        constdof[nodedof * no - 1, 0] = nodedof * no
                elif int(constrains[ii, 1]) == 1:
                    if constrains[ii, 2] == 0.0:
                        fixedof[nodedof * no - 3, 0] = nodedof * no - 2
                    else:
                        constdof[nodedof * no - 3, 0] = nodedof * no - 2
                elif int(constrains[ii, 1]) == 2:
                    if constrains[ii, 2] == 0.0:
                        fixedof[nodedof * no - 2, 0] = nodedof * no - 1
                    else:
                        constdof[nodedof * no - 2, 0] = nodedof * no - 1
                elif int(constrains[ii, 1]) == 3:
                    if constrains[ii, 2] == 0.0:
                        fixedof[nodedof * no - 1, 0] = nodedof * no
                    else:
                        constdof[nodedof * no - 1, 0] = nodedof * no

        fixedof = fixedof[np.nonzero(fixedof)]
        fixedof = fixedof - np.ones_like(fixedof)
        constdof = constdof[np.nonzero(constdof)]
        constdof = constdof - np.ones_like(constdof)
        alldof = np.arange(0, nodedof * nodetot, 1, int)
        freedof = np.setdiff1d(alldof, fixedof)
        freedof = np.setdiff1d(freedof, constdof)
        return freedof, fixedof, constdof

    # Dirichlet Non-Homogeneous
    def getDirichletNH(constrains, nodetot, nodedof):
        """
        getLoadAssembler Assembler module <ConcreteClassService>

        Returns:
            vector sparse
        """

        steps = len(np.unique(constrains[:, 3][constrains[:, 3] != 0]))
        if steps == 0:
            steps = 1
        else:
            pass
        Uc = np.zeros(
            (nodedof * nodetot, steps), dtype=np.float64
        )  # solution constrains

        for cstep in range(len(np.unique(constrains[:, 3][constrains[:, 3] != 0]))):
            forceaply = constrains[np.where(constrains[:, 3] == cstep + 1), :][0]
            nload = forceaply.shape[0]

            if nodedof == 1:
                for ii in range(nload):
                    if int(forceaply[ii, 1]) == 1:
                        gdlload = int(nodedof * forceaply[ii, 0] - (nodedof - 1))
                        Uc[gdlload, cstep] = forceaply[ii, 2]

            elif nodedof == 2:
                for ii in range(nload):
                    if int(forceaply[ii, 1]) == 1:
                        gdlload = int(nodedof * forceaply[ii, 0] - (nodedof))
                        Uc[gdlload, cstep] = forceaply[ii, 2]

                    elif int(forceaply[ii, 1]) == 2:
                        gdlload = int(nodedof * forceaply[ii, 0] - (nodedof - 1))
                        Uc[gdlload, cstep] = forceaply[ii, 2]

            elif nodedof == 3:
                for ii in range(nload):
                    if int(forceaply[ii, 1]) == 1:
                        gdlload = int(nodedof * forceaply[ii, 0] - (nodedof))
                        Uc[gdlload, cstep] = forceaply[ii, 2]

                    elif int(forceaply[ii, 1]) == 2:
                        gdlload = int(nodedof * forceaply[ii, 0] - (nodedof - 1))
                        Uc[gdlload, cstep] = forceaply[ii, 2]

                    elif int(forceaply[ii, 1]) == 3:
                        gdlload = int(nodedof * forceaply[ii, 0] - (nodedof - 2))
                        Uc[gdlload, cstep] = forceaply[ii, 2]

            else:
                pass
        return Uc
