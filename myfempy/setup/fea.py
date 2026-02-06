from __future__ import annotations

import logging
import sys
from time import time

import numpy as np
import numpy.typing as npt

from myfempy.core.utilities import setSteps
# from myfempy.core.solver import getSolver
from myfempy.io.controllers import (setElement, setGeometry, setMaterial,
                                    setMesh, setShape, setDomain, setCoupling,
                                    setPoints2NumericalIntegration)

from myfempy.plots.prevplot import preview_plot
from myfempy.setup.model import SetModel
from myfempy.setup.physics import SetPhysics
from myfempy.setup.results import setPostProcess
from myfempy.utils.utils import (clear_console, get_logo, get_version,
                                 loading_bar_v1, newDir, print_console, get_about)

__docformat__ = "google"

__doc__ = """
module for performing finite element analysis with the myfempy package.

![](docs/assets/logo2.png)
"""                             

class newAnalysis:
    """
    Setup the New Analysis to FEA simulation
    """
    def __init__(self, FEASolver: object, path: str = None):
        """Initialize a Finite Element Analysis object.

        Arguments:
            FEASolver -- analysis module to solve the problem

        Example:
            ```python
                from myfempy import newAnalysis
                from myfempy import SteadyStateLinear
                FEA = newAnalysis(SteadyStateLinear)
            ```

        Keyword Arguments:
            path -- path to output files (default: {None})
        """
        self.solver = FEASolver
        try:
            self.path = newDir(path)
        except:
            self.path = newDir("out")
        logging.basicConfig(
            filename=str(self.path) + "/" + "myfempy_api-log.log",
            encoding="utf-8",
            level=logging.DEBUG,
            filemode="w",
        )

    def Model(self, modeldata: dict) -> None:
        """
        finite element model set

        Arguments:
            modeldata -- myfempy data information

        Example:
            ```python

                mat = {
                    'NAME': 'material',
                    'EXX': 1000,
                    'VXX': 0.3
                }

                geo = {
                    'NAME':'geo',
                    'THICKN': 0.1,
                }

                nodes = [
                    [1,     0.00,   0.00,   0.00],
                    [2,     2.00,   0.00,   0.00],
                    [3,     2.00,   3.00,   0.00],
                    [4,     0.00,   2.00,   0.00],
                    [5,     0.40,   0.40,   0.00],
                    [6,     1.40,   0.60,   0.00],
                    [7,     1.50,   2.00,   0.00],
                    [8,     0.30,   1.60,   0.00],    
                        ]

                conec = [
                    [1, 1, 1, 1, 2, 6, 5],
                    [2, 1, 1, 2, 3, 7, 6],
                    [3, 1, 1, 3, 4, 8, 7],
                    [4, 1, 1, 1, 5, 8, 4],
                    [5, 1, 1, 5, 6, 7, 8],
                        ]

                modeldata = {
                    'MESH':{
                        'TYPE':'manual',
                        'COORD':nodes,
                        'INCI':conec
                    },

                    'ELEMENT':{
                        'TYPE':'structplane',
                        'SHAPE':'quad4',
                        'INTGAUSS':2,
                    },

                    'MATERIAL':{
                        'MAT':'planestress',
                        'TYPE':'isotropic',
                        'PROPMAT':[mat]
                    },

                    'GEOMETRY':{
                        'GEO':'thickness',
                        'PROPGEO':[geo]
                    }
                }

                fea.Model(modeldata)
            ```
        """
        clear_console()
        get_logo()
        print_console("mesh")
        try:
            modeldata["MESH"]["user_path"] = self.path
            Mesh = newAnalysis.__setMesh(modeldata)
            logging.info("TRY SET MESH -- SUCCESS")
        except:
            logging.warning("TRY SET MESH -- FAULT")
        try:
            Element = newAnalysis.__setElement(modeldata)
            logging.info("TRY SET ELEMENT -- SUCCESS")
        except:
            logging.warning("TRY SET ELEMENT -- FAULT")
        try:
            Shape = newAnalysis.__setShape(modeldata)
            logging.info("TRY SET SHAPE -- SUCCESS")
        except:
            logging.warning("TRY SET SHAPE -- FAULT")    
        try:
            Material = newAnalysis.__setMaterial(modeldata)
            logging.info("TRY SET MATERIAL -- SUCCESS")
        except:
            logging.warning("TRY SET MATERIAL -- FAULT")
        try:
            Geometry = newAnalysis.__setGeometry(modeldata)
            logging.info("TRY SET GEOMETRY -- SUCCESS")
        except:
            logging.warning("TRY SET GEOMETRY -- FAULT")
        try:
            GaussPoints = newAnalysis.__setIntGauss(modeldata)
            logging.info("TRY SET GAUSS POINTS -- SUCCESS")
        except:
            logging.warning("TRY SET GAUSS POINTS -- FAULT")            
        try:
            self.model = SetModel(Mesh, Element, Shape, Material, Geometry)
            self.model.modeldata = modeldata
            # self.model.intgauss = GaussPoints
            logging.info("TRY SET FEMODEL -- SUCCESS")
        except:
            logging.warning("TRY SET FEMODEL -- FAULT")
            
        
        self.model.inci = newAnalysis.getInci(self)
        self.model.coord = newAnalysis.getCoord(self)
        self.model.tabmat = newAnalysis.getTabmat(self)
        self.model.tabgeo = newAnalysis.getTabgeo(self)
        self.model.intgauss = GaussPoints

        self.model.modelinfo = dict()
        try:
            self.model.regions = newAnalysis.getRegions(self)
        except:
            self.model.regions = []
        elem_set = self.model.element.getElementSet()
        self.model.modelinfo["tensor"] = len(elem_set["tensor"])
        self.model.modelinfo["dofs"] = elem_set["dofs"]
        self.model.modelinfo["nodedof"] = len(elem_set["dofs"]["d"])
        self.model.modelinfo["type_element"] = elem_set["key"]
        shape_set = self.model.shape.getShapeSet()
        self.model.modelinfo["shapeid"] = shape_set["id"]
        self.model.modelinfo["nodecon"] = len(shape_set["nodes"])
        self.model.modelinfo["elemdofs"] = len(shape_set["nodes"]) * self.model.modelinfo["nodedof"]
        self.model.modelinfo["type_shape"] = shape_set["key"]
        self.model.modelinfo["elemid"] = int(f'{elem_set["id"]}{shape_set["id"]}')
        self.model.modelinfo["nnode"] = len(self.model.coord)
        self.model.modelinfo["nelem"] = len(self.model.inci)
        self.model.modelinfo["fulldofs"] = len(elem_set["dofs"]["d"]) * len(self.model.coord)

        self.model.elemvol = newAnalysis.getElementVolume(
            self,
            self.model.inci,
            self.model.coord,
            self.model.tabgeo,
        )

    def Physic(self, physicdata: dict) -> None:
        """
        set load and boundary conditions

        Arguments:
            physicdata -- myfempy data information

        Example:
            ```python
                bcfix = {
                    'TYPE':'fixed',
                    'DOF':'full',
                    'DIR':'node',
                    'LOC':{'x':0, 'y':0, 'z':0},
                }

                def displ_patch_test(x,y):
                    ux = 0.002*x
                    uy = -0.0006*y
                    return ux, uy

                bcnh_node1_ux = {
                    'TYPE':'displ',
                    'DOF':'ux',
                    'DIR':'node',
                    'LOC':{'x':0.0, 'y':0, 'z':0},
                    'VAL':[displ_patch_test(0.0,0.00)[0]]
                } 

                bcnh_node1_uy = {
                    'TYPE':'displ',
                    'DOF':'uy',
                    'DIR':'node',
                    'LOC':{'x':0.0, 'y':0, 'z':0},
                    'VAL':[displ_patch_test(0.0,0.00)[1]]
                } 

                bcnh_node2_ux = {
                    'TYPE':'displ',
                    'DOF':'ux',
                    'DIR':'node',
                    'LOC':{'x':2.0, 'y':0, 'z':0},
                    'VAL':[displ_patch_test(2.0,0.00)[0]]
                } 

                bcnh_node2_uy = {
                    'TYPE':'displ',
                    'DOF':'uy',
                    'DIR':'node',
                    'LOC':{'x':2.0, 'y':0, 'z':0},
                    'VAL':[displ_patch_test(2.0,0.00)[1]]
                } 

                bcnh_node3_ux = {
                    'TYPE':'displ',
                    'DOF':'ux',
                    'DIR':'node',
                    'LOC':{'x':2.0, 'y':3.0,'z':0},
                    'VAL':[displ_patch_test(2.0,3.0)[0]]
                } 

                bcnh_node3_uy = {
                    'TYPE':'displ',
                    'DOF':'uy',
                    'DIR':'node',
                    'LOC':{'x':2.0, 'y':3.0, 'z':0},
                    'VAL':[displ_patch_test(2.0,3.0)[1]]
                } 


                bcnh_node4_ux = {
                    'TYPE':'displ',
                    'DOF':'ux',
                    'DIR':'node',
                    'LOC':{'x':0, 'y':2.0,'z':0},
                    'VAL':[displ_patch_test(0,2.0)[0]]
                } 

                bcnh_node4_uy = {
                    'TYPE':'displ',
                    'DOF':'uy',
                    'DIR':'node',
                    'LOC':{'x':0, 'y':2.0, 'z':0},
                    'VAL':[displ_patch_test(0,2.0)[1]]
                } 

                physicdata = {
                    'PHYSIC':{
                        'DOMAIN':'structural',
                        'LOAD':[],
                        'BOUNDCOND':[
                                    # bcfix,
                                    bcnh_node1_ux, bcnh_node1_uy,
                                    bcnh_node2_ux, bcnh_node2_uy,
                                    bcnh_node3_ux, bcnh_node3_uy,
                                    bcnh_node4_ux, bcnh_node4_uy,]
                    }
                }
                fea.Physic(physicdata)

            ```
        """
        print_console("pre")
        # self.model.modelinfo["physic"] = physicdata["PHYSIC"]
        try:
            Loads, BoundCond = newAnalysis.__setDomain(physicdata)
            self.physic = SetPhysics(self.model, Loads, BoundCond)
            self.physic.physicdata = physicdata
            logging.info("TRY SET PHYSICS -- SUCCESS")
        except:
            self.physic = []
            logging.warning("TRY SET PHYSICS -- FAULT")

        try:
            self.physic.forces = newAnalysis.getLoadApply(self)
            logging.info("TRY SET PHYSICS.FORCES -- SUCCESS")
        except:
            self.physic.forces = []
            logging.warning("TRY SET PHYSICS.FORCES -- FAULT")

        if "COUPLING" in physicdata.keys():
            coupling_load_zero = self.physic.forces
            LoadCoup, BoundCond = newAnalysis.__setCoupling(physicdata)
            self.physic = SetPhysics(self.model, LoadCoup, BoundCond)
            self.physic.physicdata = physicdata
            self.physic.forces = np.append(
                coupling_load_zero, newAnalysis.getCouplingInterface(self), axis=0
            )
            
        try:
            constrains = newAnalysis.getBCApply(self)
            if any(constrains[:, 1] == 11):
                self.physic.csleft = constrains[
                    np.where(constrains[:, 1] == 11)[0], 0
                ]
                self.physic.csright = constrains[
                    np.where(constrains[:, 1] == 12)[0], 0
                ]
                self.physic.constrains = constrains
            else:
                self.physic.constrains = constrains
            logging.info("TRY SET PHYSICS.CONSTRAINS -- SUCCESS")
        except:
            self.physic.constrains = []
            logging.warning("TRY SET PHYSICS.CONSTRAINS -- FAULT")

        # self.loadaply = FEANewAnalysis.getLoadApply(self)
        # self.constrains = FEANewAnalysis.getBCApply(self)

    def Assembly(self, Model: object) -> npt.NDArray[np.float64]:
        """assembly of fe model algebric system

        Arguments:
            Model -- analysis model object

        Returns:
            matrix -- npt.NDArray[np.float64]
            forcelist --npt.NDArray[np.float64]
        """
        try:
            matrix = newAnalysis.getGlobalMatrix(self, Model, self.model.inci, self.model.coord, self.model.tabmat, self.model.tabgeo, self.model.intgauss, self.symm, self.mp)
            logging.info("TRY RUN GLOBAL ASSEMBLY -- SUCCESS")
        except:
            logging.warning("TRY RUN GLOBAL ASSEMBLY -- FAULT")
        loadaply = self.physic.forces
        try:
            matrix = newAnalysis.getUpdateMatrix(self, matrix, loadaply)
            logging.info("TRY RUN UPDATE ASSEMBLY -- SUCCESS")
        except:
            logging.warning("TRY RUN UPDATE ASSEMBLY -- FAULT")
        try:
            forcelist = newAnalysis.getLoadArray(self, loadaply)
            logging.info("TRY RUN LOAD ASSEMBLY -- SUCCESS")
        except:
            logging.warning("TRY RUN LOAD ASSEMBLY -- FAULT")
        
        return matrix, forcelist
    # 
    def Solve(self, solverset=None) -> dict:
        """run the solver set

        Keyword Arguments:
            solverset -- configuration file for solution (default: {None})

        Example:
            ```python
                solverset = {
                    'STEPSET':{
                        'type':'table',
                        'start':0,
                        'end':1,
                        'step':1
                    },
                    'SYMM':False,
                    'MP':False
                }

                solverdata = fea.Solve(solverset)
            ```

        Returns:
            solverset
        """
        print_console("solver")
        try:
            solverset["solverstatus"] = dict()
            self.symm = solverset["SYMM"]
            if self.symm:
                solverset["solverstatus"]["typeasmb"] = "SYMMETRIC"
            else:
                solverset["solverstatus"]["typeasmb"] = "FULL"
        except:
            self.symm = False
            solverset["solverstatus"]["typeasmb"] = "FULL"
        try:
            self.mp = solverset["MP"]
            solverset["solverstatus"]["ncpu"] = (
                "PARALLEL_" + str(solverset["MP"]) + "_CORES"
            )
        except:
            self.mp = 0
            solverset["solverstatus"]["ncpu"] = "SERIAL_" + str(1) + "_CORE"
        # loading_bar_v1(10,"SOLVER")
        starttime = time()
        assembly, forcelist = newAnalysis.Assembly(self, Model=self.model)
        endttime = time()
        solverset["solverstatus"]["timeasb"] = abs(endttime - starttime)
        solverset["solverstatus"]["memorysize"] = (assembly["stiffness"].todense().nbytes)/1e6
        # loading_bar_v1(50,"SOLVER")
        try:
            constrains = self.physic.constrains
            freedof, fixedof, constdof = newAnalysis.getConstrains(self, constrains)
            logging.info("TRY RUN CONSTRAINS -- SUCCESS")
        except:
            logging.warning("TRY RUN CONSTRAINS -- FAULT")
        # loading_bar_v1(60,"SOLVER")
        constrainsdof = dict()
        constrainsdof["freedof"] = freedof
        constrainsdof["fixedof"] = fixedof
        constrainsdof["constdof"] = constdof
        try:
            Uc = newAnalysis.getDirichletNH(self, constrains)
            logging.info("TRY SET DNH CONSTRAINS -- SUCCESS")
        except:
            logging.warning("TRY SET DNH CONSTRAINS -- FAULT")
        nsteps = setSteps(solverset["STEPSET"])
        if forcelist.shape[1] != nsteps:
            forcelist = np.repeat(forcelist, nsteps, axis=1)
        else:
            pass
        assembly["loads"] = forcelist
        if Uc.shape[1] != nsteps:
            Uc = np.repeat(Uc, nsteps, axis=1)
        else:
            pass
        assembly["bcdirnh"] = Uc
        # loading_bar_v1(80,"SOLVER")
        try:
            starttime = time()
            solverset["solution"] = self.solver.runSolve(self.model, self.physic, assembly, constrainsdof, solverset)
            endttime = time()
            solverset["solverstatus"]["timesim"] = abs(endttime - starttime)
            logging.info("TRY RUN SOLVER -- SUCCESS")
        except:
            logging.warning("TRY RUN SOLVER -- FAULT")
        # loading_bar_v1(100,"SOLVER")
        return solverset

    def PreviewAnalysis(self, previewdata) -> None:
        """
        preview the model+physic set

        Arguments:
            previewdata -- myfempy data information

        Example:
            ```python
                previewset = {
                    'RENDER':{
                        'filename':'patchtest',
                        'show':True,
                        'scale':5,
                        'savepng':True,
                        'lines':True
                    }
                }
                fea.PreviewAnalysis(previewset)
            ```

        """        
        try:
            preview_plot(self.model, previewdata, str(self.path), self.physic)
            logging.info("TRY RUN PREVIEW PLOT -- SUCCESS")
        except:
            preview_plot(self.model, previewdata, str(self.path))
            logging.warning("TRY RUN PREVIEW PLOT -- FAULT")

    def PostProcess(self, postprocset) -> dict:
        """
        prost process the solution

        Arguments:
            postprocset -- configuration file

        Example:
            ```python
                postprocset = {"SOLVERDATA": solverdata,
                                "COMPUTER": {'structural': {'displ': True, 'stress': True}},
                                "PLOTSET": {'show': True, 'filename': 'PatchTest', 'savepng': True},
                                "OUTPUT": {'log': True, 'get':{
                                        'nelem': True,
                                        'nnode': True,
                                        'inci': True,
                                        'coord':True,
                                        'tabmat':True,
                                        'tabgeo':True,
                                        'boundcond_list':True,
                                        'forces_list':True,
                                    }
                            }}
                postprocdata = fea.PostProcess(postprocset)
            ```

        Returns:
            postprocdata
        """
        print_console("post")
        postprocdata = []
        try:
            if "COMPUTER" in postprocset.keys():
                postprocdata = setPostProcess.getCompute(self, postprocset)
            if "TRACKER" in postprocset.keys():
                setPostProcess.getTracker(self, postprocset, postprocdata)
            if "REPORT" in postprocset.keys():
                postprocdata["log"] = []
                log_file = setPostProcess.getLog(self, postprocset, postprocdata)
                postprocdata["log"].append(log_file)
            logging.info("TRY GET POST PROCESS -- SUCCESS")
        except:
                logging.warning("TRY GET POST PROCESS -- FAULT")
        return postprocdata

    # GET MODEL
    def getModel(self) -> object:
        """get object related to the current simulation model

        Returns:
            object
        """
        return self.model

    def getModelInfo(self) -> dict:
        """get model info

        Returns:
            dict
        """
        return self.model.modelinfo

    def getInci(self) -> npt.NDArray[np.float64]:
        """get mesh properties table

        Returns:
            npt.NDArray[np.float64]
        """
        return self.model.getInci(self.model.modeldata)

    def getCoord(self) -> npt.NDArray[np.float64]:
        """get mesh grid coordinate

        Returns:
            npt.NDArray[np.float64]
        """
        return self.model.getCoord(self.model.modeldata)

    def getTabmat(self) -> list:
        """get table of material properties

        Returns:
            list
        """
        return self.model.getTabMat(self.model.modeldata)

    def getTabgeo(self) -> list:
        """get table of geometry properties

        Returns:
            list
        """
        return self.model.getTabGeo(self.model.modeldata)

    def getIntGauss(self) -> int:
        """get Gaussian numerical integration number

        Returns:
            int
        """
        return self.model.intgauss

    def getElementVolume(self, inci:npt.NDArray[np.float64], coord:npt.NDArray[np.float64], tabgeo:list) -> npt.NDArray[np.float64]:
        """get elements volumes list

        Arguments:
            inci -- npt.NDArray[np.float64]
            coord -- npt.NDArray[np.float64]
            tabgeo -- list

        Returns:
            npt.NDArray[np.float64]
        """
        vol = np.zeros((inci.shape[0]))
        for ee in range(inci.shape[0]):
            vol[ee] = self.model.element.getElementVolume(
                self.model, inci, coord, tabgeo, ee
            )
        return vol

    def getElemStifLinearMat(
        self, inci: npt.NDArray[np.float64], coord: npt.NDArray[np.float64], tabmat: list, tabgeo: list, intgauss: int, element_number: int
    ) -> npt.NDArray[np.float64]:
        """get elementar linear stiffness matrix

        Arguments:
            inci -- npt.NDArray[np.float64]
            coord -- npt.NDArray[np.float64]
            tabmat -- list
            tabgeo -- list
            intgauss -- int
            element_number -- int

        Returns:
            npt.NDArray[np.float64]
        """
        return self.model.element.getStifLinearMat(
            self.model, inci, coord, tabmat, tabgeo, intgauss, element_number
        )

    def getElemMassConsistentMat(
        self, inci: npt.NDArray[np.float64], coord: npt.NDArray[np.float64], tabmat: list, tabgeo: list, intgauss: int, element_number: int
    ) -> npt.NDArray[np.float64]:
        """get elementar linear mass matrix

        Arguments:
            inci -- npt.NDArray[np.float64]
            coord -- npt.NDArray[np.float64]
            tabmat -- list
            tabgeo -- list
            intgauss -- int
            element_number -- int

        Returns:
            npt.NDArray[np.float64]
        """
        return self.model.element.getMassConsistentMat(
            self.model, inci, coord, tabmat, tabgeo, intgauss, element_number
        )
    
    def getRegions(self) -> list:
        """get regions from gmsh mesh only

        Returns:
            list
        """
        return self.model.mesh.getRegionsList(
            self.model.mesh.getElementConection(self.model.modeldata["MESH"])
        )

    # GET SOLVER
    def getGlobalMatrix(self, Model, inci:npt.NDArray[np.float64] = None, coord:npt.NDArray[np.float64] = None, tabmat:list = None, tabgeo:list = None, intgauss:int = None, SYMM:bool=None, MP:bool=None) -> npt.NDArray[np.float64]:
        """get global assembler matrix 

        Arguments:
            Model -- object related to the current simulation model

        Keyword Arguments:
            inci -- mesh properties table (default: {None})
            coord -- grid coordinate table (default: {None})
            tabmat -- table of material properties (default: {None})
            tabgeo -- table of geometry properties (default: {None})
            intgauss -- Gaussian numerical integration number (default: {None})
            SYMM -- symmetric assembler (default: {None})
            MP -- multi-processing (default: {None})

        Returns:
            npt.NDArray[np.float64]
        """
        return self.solver.getMatrixAssembler(Model, inci = inci, coord = coord, tabmat = tabmat, tabgeo = tabgeo, intgauss = intgauss, SYMM=SYMM, MP=MP)

    def getConstrains(self, constrains:list) -> npt.NDArray[np.float64]:
        """get constrains boundary conditions list

        Arguments:
            constrains -- _description_

        Returns:
            _description_
        """
        nodetot = len(self.model.coord)
        return self.solver.getConstrains(
            constrains, nodetot, self.model.modelinfo["nodedof"]
        )

    def getDirichletNH(self, constrains:list) -> npt.NDArray[np.float64]:
        """get dirichlet non-homogeneous boundary conditions list

        Arguments:
            constrains -- _description_

        Returns:
            _description_
        """
        nodetot = len(self.model.coord)
        return self.solver.getDirichletNH(
            constrains, nodetot, self.model.modelinfo["nodedof"]
        )

    def getLoadArray(self, loadaply:list) -> npt.NDArray[np.float64]:
        """get loads arrays

        Arguments:
            loadaply -- _description_

        Returns:
            _description_
        """
        nodetot = len(self.model.coord)
        return self.solver.getLoadAssembler(
            loadaply, nodetot, self.model.modelinfo["nodedof"]
        )

    # GET PHYSIC
    def getPhysic(self) -> object:
        """get physic object

        Returns:
            _description_
        """
        return self.physic

    def getForceList(self) -> list:
        """get forces list

        Returns:
            _description_
        """
        return self.physic.getForceList(self.modelinfo["domain"])

    def getBoundCondList(self) -> list:
        """get boundary conditions list

        Returns:
            _description_
        """
        return self.physic.getBoundCondList(self.modelinfo["domain"])

    def getLoadApply(self) -> npt.NDArray[np.float64]:
        """get loads applied array on model

        Returns:
            _description_
        """
        return self.physic.getLoadApply(self.physic.physicdata)

    def getBCApply(self) -> npt.NDArray[np.float64]:
        """get boundary conditions applied array on model

        Returns:
            _description_
        """
        return self.physic.getBoundCondApply(self.physic.physicdata)
    
    def getCouplingInterface(self) -> list:
        """get coupling interface list

        Returns:
            _description_
        """
        return self.physic.getLoadCoup(self.physic.physicdata)

    def getUpdateMatrix(self, matrix, addval) -> npt.NDArray[np.float64]:
        """get updated matrix

        Arguments:
            matrix -- _description_
            addval -- _description_

        Returns:
            _description_
        """
        return self.physic.getUpdateMatrix(matrix, addval)
    
    def getElementFromNodesList(self, nodelist) -> list:
        """get elements from nodes list

        Arguments:
            nodelist -- _description_

        Returns:
            _description_
        """
        return self.physic.getElementList(self.model.inci, nodelist)

    # OUTHERS
    def getNodesFromRegions(self, set: int, type: str) -> list:
        """get nodes from regions gmsh mesh only

        Arguments:
            set -- _description_
            type -- _description_

        Returns:
            _description_
        """
        if type == "point":
            domain_nodelist = newAnalysis.getRegions(self)[0][1][set - 1][1]
        elif type == "line":
            domain_nodelist = newAnalysis.getRegions(self)[1][1][set - 1][1]
        elif type == "plane":
            domain_nodelist = newAnalysis.getRegions(self)[2][1][set - 1][1]
        else:
            domain_nodelist = []
        return domain_nodelist
    
    def getAbout(self):
        get_about()

    # settings FEA ANALYSIS <privates>
    def __setMesh(modeldata):
        set_mesh = dict()
        set_mesh = modeldata["MESH"]
        set_mesh["SHAPE"] = modeldata["ELEMENT"]["SHAPE"]
        set_mesh["user_path"] = modeldata["MESH"]["user_path"]
        return setMesh(set_mesh)

    def __setShape(modeldata):
        return setShape(modeldata["ELEMENT"])

    def __setElement(modeldata):
        return setElement(modeldata["ELEMENT"])

    def __setMaterial(modeldata):
        return setMaterial(modeldata["MATERIAL"])

    def __setGeometry(modeldata):
        return setGeometry(modeldata["GEOMETRY"])

    def __setDomain(modeldata):
        return setDomain(modeldata["PHYSIC"])

    def __setCoupling(modeldata):
        return setCoupling(modeldata["COUPLING"])
    
    def __setIntGauss(modeldata):
        if "INTGAUSS" in modeldata["ELEMENT"].keys():
            intgauss = modeldata["ELEMENT"]["INTGAUSS"]
        else:
            intgauss = setPoints2NumericalIntegration(modeldata["ELEMENT"]["SHAPE"])
        return intgauss


    # def __setSolution(self, solvedata):
    #     # self.inci = FEANewAnalysis.getInci(self)
    #     # self.coord = FEANewAnalysis.getCoord(self)
    #     # self.tabmat = FEANewAnalysis.getTabmat(self)
    #     # self.tabgeo = FEANewAnalysis.getTabgeo(self)
    #     # self.intgauss = FEANewAnalysis.getIntGauss(self)
    #     # elem_set = self.model.element.getElementSet()
    #     # nodedof = len(elem_set["dofs"])
    #     self.solver.fulldofs = (self.modelinfo['dofs']) * len(self.modelinfo['nnode'])
    #     self.solver.solvedata = solvedata
