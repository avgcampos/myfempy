from __future__ import annotations


from numpy import dot, float64, zeros, array, ix_, where, arange, concatenate, in1d
from scipy.sparse.linalg import spsolve

from myfempy.core.solver.assemblerfull import AssemblerFULL
from myfempy.core.solver.assemblerfull_parallel import AssemblerFULLPOOL
from myfempy.core.solver.assemblersymm import AssemblerSYMM
# from myfempy.core.alglin import linsolve_spsolve
from myfempy.core.solver.solver import Solver
from myfempy.core.utilities import setSteps, gauss_points

class HomogenPlane(Solver):
    """
    Homogenization Plane Solver Class <ConcreteClassService>
    """

    # @profile
    def getMatrixAssembler(Model, inci = None, coord = None, tabmat = None, tabgeo = None, intgauss = None, SYMM=None, MP=None):
        matrix = dict()
        if SYMM:
            matrix["stiffness"] = AssemblerSYMM.getLinearStiffnessGlobalMatrixAssembler(Model, inci, coord, tabmat, tabgeo, intgauss)
        else:
            if MP:
                matrix["stiffness"] = AssemblerFULLPOOL.getLinearStiffnessGlobalMatrixAssembler(Model, inci, coord, tabmat, tabgeo, intgauss, MP)
            else:
                matrix["stiffness"] = AssemblerFULL.getLinearStiffnessGlobalMatrixAssembler(Model, inci, coord, tabmat, tabgeo, intgauss)
        return matrix

    def getLoadAssembler(loadaply, nodetot, nodedof):
        return AssemblerFULL.getLoadAssembler(loadaply, nodetot, nodedof)

    def getConstrains(constrains, nodetot, nodedof):
        nodes_constrain_XX = constrains[where(constrains[:, 3] == 1)]
        nodes_constrain_YY = constrains[where(constrains[:, 3] == 2)]
        nodes_constrain_XY = constrains[where(constrains[:, 3] == 3)]

        __, dofs_constrain_XX, __ = AssemblerFULL.getConstrains(nodes_constrain_XX, nodetot, nodedof)
        __, dofs_constrain_YY, __ = AssemblerFULL.getConstrains(nodes_constrain_YY, nodetot, nodedof)
        __, dofs_constrain_XY, __ = AssemblerFULL.getConstrains(nodes_constrain_XY, nodetot, nodedof)

        freedof = []
        fixedof = [dofs_constrain_XX, dofs_constrain_YY, dofs_constrain_XY]
        constdof = []

        return freedof, fixedof, constdof

    def getDirichletNH(constrains, nodetot, nodedof):
        return AssemblerFULL.getDirichletNH(constrains, nodetot, nodedof)

    def runSolve(Model, Physic, assembly, constrainsdof, solverset):
        solution = dict()
        elem_set = Model.element.getElementSet()
        nodedof = len(elem_set["dofs"]["d"])
        ntensor = len(elem_set['tensor'])
        shape_set = Model.shape.getShapeSet()
        type_shape = shape_set["key"]
        ndofs = Model.modelinfo["fulldofs"]

        fulldofs = arange(0, ndofs, 1, int)

        stiffness = assembly["stiffness"]
        forcelist = assembly["loads"]

        fixeddof_case_epsXX = constrainsdof["fixedof"][0]
        fixeddof_case_epsYY = constrainsdof["fixedof"][1]
        fixeddof_case_epsXY = constrainsdof["fixedof"][2]

        freedof_pc_case_epsXX = where(in1d(fulldofs, fixeddof_case_epsXX, assume_unique=True) == False)[0]
        freedof_pc_case_epsYY = where(in1d(fulldofs, fixeddof_case_epsYY, assume_unique=True) == False)[0]
        freedof_pc_case_epsXY = where(in1d(fulldofs, fixeddof_case_epsXY, assume_unique=True) == False)[0]

        freedof_pc = [freedof_pc_case_epsXX, freedof_pc_case_epsYY, freedof_pc_case_epsXY]

        U = zeros((ndofs, ntensor), dtype=float64)    # empty((fulldofs, nsteps))
        for step in range(ntensor):
            try:
                X = spsolve(stiffness[:, freedof_pc[step]][freedof_pc[step], :], forcelist[freedof_pc[step], step])
            except:
                pass
            U[freedof_pc[step], step] = X
        
        inci = Model.inci
        coord = Model.coord
        tabmat = Model.tabmat
        tabgeo = Model.tabgeo
        intgauss = Model.intgauss

        CH = zeros((ntensor, ntensor))
        for elm in range(inci.shape[0]):
            nodelist = Model.shape.getNodeList(inci, elm)
            elementcoord = Model.shape.getNodeCoord(coord, nodelist)
            t = tabgeo[int(inci[elm, 3] - 1)]["THICKN"]
            Ci = Model.material.getElasticTensor(tabmat, inci, elm)
            pt, wt = gauss_points(type_shape, intgauss)
            loc = Model.shape.getLocKey(nodelist, nodedof)
            ui = U[ix_(loc)]
            CHelm = zeros((ntensor, ntensor))
            for ip in range(intgauss):
                for jp in range(intgauss):
                    detJ = Model.shape.getdetJacobi(array([pt[ip], pt[jp]]), elementcoord)
                    diffN = Model.shape.getDiffShapeFuntion(array([pt[ip], pt[jp]]), nodedof)
                    invJ = Model.shape.getinvJacobi(array([pt[ip], pt[jp]]), elementcoord, nodedof)
                    B = Model.element.getB(diffN, invJ)
                    CHelm +=  (Ci - dot(Ci, dot(B, ui))) * t * abs(detJ) * wt[ip] * wt[jp]

            CH += CHelm

        Yx = max(coord[:,1])
        Yy = max(coord[:,2])

        CH = CH/(Yx * Yy * t)

        solution["U"] = U
        solution["CH"] = CH

        return solution
